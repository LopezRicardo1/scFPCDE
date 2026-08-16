#' Run the scFPCDE Analysis Pipeline
#'
#' Runs the complete workflow: ordering observations by pseudotime, centering or
#' scaling expression, estimating an initial FPCA basis, refining the basis with
#' the strongest gene trajectories, and applying permutation tests.
#'
#' @param yt Numeric expression matrix with cells in rows and genes in columns.
#'   Gene names are required as column names.
#' @param tt Numeric pseudotime vector with one value per row of `yt`.
#' @param L Positive integer giving the number of harmonics to retain.
#' @param r_pen Non-negative smoothing penalty passed to [fda::fdPar()].
#' @param nbasis Integer number of B-spline basis functions; must be at least 4.
#' @param n_perm Positive integer giving the number of permutations per test.
#' @param topvarper Proportion in `(0, 1]` of genes used to refine the FPCA
#'   eigenfunctions.
#' @param center Logical; center each gene before fitting.
#' @param scale Logical; scale each gene to unit variance before fitting.
#' @param ncores Positive integer giving the number of parallel workers. Use
#'   `ncores = 1` for sequential execution.
#' @param use_FPC_F Logical; if `TRUE`, also run [scFPCDE_F_test()].
#' @param fpc_varmax Logical; if `TRUE`, apply VARIMAX rotation to FPCA
#'   harmonics in both fitting steps.
#'
#' @return A list with components:
#' \describe{
#'   \item{fpca_result}{The refined result from [scFPCDE_fit_fpca()].}
#'   \item{D_test_result}{The data frame returned by [scFPCDE_D_test()].}
#'   \item{F_test_result}{The data frame returned by [scFPCDE_F_test()] when
#'   `use_FPC_F = TRUE`; otherwise `NULL`.}
#' }
#'
#' @details
#' The fitted values in `fpca_result$xt_hat` correspond to the centered or
#' scaled matrix used internally. For deterministic permutations, call
#' [set.seed()] and use `ncores = 1`.
#'
#' @examples
#' \donttest{
#' data(scFPCDE_simdata)
#' cells <- seq(1, 1000, length.out = 100)
#' cells <- unique(round(cells))
#' genes <- c(1:20, 101:120)
#' set.seed(1)
#' result <- scFPCDE_run(
#'   scFPCDE_simdata$yt[cells, genes],
#'   scFPCDE_simdata$tt[cells],
#'   nbasis = 20,
#'   n_perm = 10,
#'   ncores = 1
#' )
#' head(result$D_test_result)
#' }
#' @export
scFPCDE_run <- function(yt, tt, L = 2, r_pen = 1e-3, nbasis = 50,
                        n_perm = 1000, topvarper = 0.1,
                        center = TRUE, scale = FALSE,
                        ncores = 2, use_FPC_F = FALSE,
                        fpc_varmax = TRUE) {

  scFPCDE_validate_expression(yt, tt)
  scFPCDE_validate_fpca_arguments(L, r_pen, nbasis, fpc_varmax)
  n_perm <- scFPCDE_validate_integer(n_perm, "n_perm")
  ncores <- scFPCDE_validate_integer(ncores, "ncores")
  scFPCDE_validate_probability(topvarper, "topvarper")
  if (length(center) != 1L || is.na(center) || !is.logical(center)) {
    stop("`center` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (length(scale) != 1L || is.na(scale) || !is.logical(scale)) {
    stop("`scale` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (length(use_FPC_F) != 1L || is.na(use_FPC_F) || !is.logical(use_FPC_F)) {
    stop("`use_FPC_F` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  # Order by pseudotime to ensure consistent trajectory alignment
  ord <- order(tt)
  tt <- tt[ord]
  yt <- yt[ord, , drop = FALSE]

  # Preprocess gene expression: center and/or scale each gene
  yt <- scale(yt, center = center, scale = scale)

  # Initial FPCA using all genes
  fpca_result_full <- scFPCDE_fit_fpca(
    yt, tt, L = L, r_pen = r_pen, nbasis = nbasis,
    fpc_varmax = fpc_varmax
  )

  # Select top variable genes
  D_stat <- rowSums(fpca_result_full$scores^2)
  n_genes <- ncol(yt)
  top_n <- ceiling(topvarper * n_genes)
  top_idx <- utils::head(order(D_stat, decreasing = TRUE), top_n)

  # Refit FPCA on top variable genes
  fpca_result_top <- scFPCDE_fit_fpca(
    yt, tt, L = L, r_pen = r_pen, nbasis = nbasis,
    topvarsub = top_idx,
    fpc_varmax = fpc_varmax
  )

  # Run D-test
  D_test_result <- scFPCDE_D_test(
    yt, tt, fpca_result_top, n_perm = n_perm, ncores = ncores
  )

  # Conditionally run F-test
  F_test_result <- NULL
  if (use_FPC_F) {
    F_test_result <- scFPCDE_F_test(
      yt, tt, fpca_result_top, n_perm = n_perm, ncores = ncores
    )
  }

  list(
    fpca_result = fpca_result_top,
    D_test_result = D_test_result,
    F_test_result = F_test_result
  )
}
