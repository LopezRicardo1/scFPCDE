# Build the packaged HB6 trajectory inputs from the exact post-processing
# objects used by the manuscript analysis.
#
# Usage:
#   Rscript data-raw/scFPCDE_hb6.R /path/to/bcell_two_traj
#
# Alternatively, set SCFPCDE_HB6_RESULTS_DIR to that directory. It must contain
# res_traj1.rds and res_traj2.rds. The source objects are read sequentially so
# that both multi-gigabyte analysis objects are never retained at the same time.

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  Sys.getenv("SCFPCDE_HB6_RESULTS_DIR", unset = "")
}

if (!nzchar(results_dir)) {
  stop(
    "Supply the bcell_two_traj directory as the first argument or set ",
    "SCFPCDE_HB6_RESULTS_DIR.",
    call. = FALSE
  )
}

source_files <- file.path(results_dir, c("res_traj1.rds", "res_traj2.rds"))
if (any(!file.exists(source_files))) {
  stop(
    "Missing required source file(s): ",
    paste(source_files[!file.exists(source_files)], collapse = ", "),
    call. = FALSE
  )
}

extract_trajectory <- function(path, trajectory_label) {
  result <- readRDS(path)
  required <- c("config", "cds_sub", "y", "tt", "tt.cluster")
  missing_fields <- setdiff(required, names(result))
  if (length(missing_fields)) {
    stop(
      basename(path), " is missing: ", paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  cell_id <- names(result$tt)
  gene_id <- unique(colnames(result$y))
  if (is.null(cell_id) || is.null(gene_id)) {
    stop("Pseudotime and analysis genes must be named.", call. = FALSE)
  }
  if (!all(cell_id %in% colnames(result$cds_sub))) {
    stop("Some paper-analysis cells are absent from cds_sub.", call. = FALSE)
  }
  if (!all(gene_id %in% rownames(result$cds_sub))) {
    stop("Some paper-analysis genes are absent from cds_sub.", call. = FALSE)
  }

  logcounts <- SummarizedExperiment::assay(result$cds_sub, "logcounts")
  yt <- t(as.matrix(logcounts[gene_id, cell_id, drop = FALSE]))
  tt <- as.numeric(result$tt[cell_id])
  names(tt) <- cell_id
  clusters <- droplevels(result$tt.cluster[cell_id])
  names(clusters) <- cell_id

  config <- result$config
  trajectory <- list(
    yt = yt,
    tt = tt,
    clusters = clusters,
    cell_id = cell_id,
    gene_id = gene_id,
    metadata = list(
      trajectory = trajectory_label,
      analysis_tag = config$analysis_tag,
      requested_clusters = config$clusters,
      source_file = basename(path),
      source_object = "cds_sub",
      expression_assay = "logcounts",
      expression_centered = FALSE,
      gene_selection = paste(
        "Unique gene columns from the official paper analysis y matrix;",
        "the duplicated ACTB column in trajectory 1 is retained once."
      ),
      paper_parameters = list(
        L = config$L,
        nbasis = config$nbasis_range,
        roughness_grid = config$r_pen_range_local,
        topvarper_tune = config$topvarper_tune,
        topvarper_run = config$topvarper_run,
        n_perm = config$n_perm_main
      ),
      study = list(
        donor = "HB6",
        accession = "E-MTAB-9544",
        doi = "10.3389/fimmu.2021.602539"
      )
    )
  )

  stopifnot(
    identical(rownames(trajectory$yt), trajectory$cell_id),
    identical(colnames(trajectory$yt), trajectory$gene_id),
    identical(names(trajectory$tt), trajectory$cell_id),
    identical(names(trajectory$clusters), trajectory$cell_id),
    length(unique(trajectory$gene_id)) == length(trajectory$gene_id),
    all(is.finite(trajectory$yt)),
    all(is.finite(trajectory$tt)),
    !anyNA(trajectory$clusters)
  )

  rm(result, logcounts)
  invisible(gc())
  trajectory
}

scFPCDE_hb6 <- list(
  traj1 = extract_trajectory(source_files[[1L]], "Trans - Naive - C-mem1"),
  traj2 = extract_trajectory(source_files[[2L]], "M-mem1 - DN2 - DN3")
)

dir.create("data", showWarnings = FALSE)
save(scFPCDE_hb6, file = "data/scFPCDE_hb6.rda", compress = "xz")

message(
  "Saved data/scFPCDE_hb6.rda (",
  paste(vapply(scFPCDE_hb6, function(x) {
    paste0(nrow(x$yt), " cells x ", ncol(x$yt), " genes")
  }, character(1)), collapse = "; "),
  ")."
)
