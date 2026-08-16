#' Tune FPCA Smoothing Parameters by Generalized Cross-Validation
#'
#' Evaluates a grid of smoothing penalties and B-spline basis sizes using mean
#' generalized cross-validation (GCV) over genes selected from an initial FPCA
#' fit.
#'
#' @param yt Numeric expression matrix with observations in rows and genes in
#'   columns.
#' @param tt Numeric pseudotime vector with one value per row of `yt`.
#' @param L Positive integer giving the number of harmonics to retain.
#' @param r_pen_range Numeric vector of non-negative smoothing penalties.
#' @param nbasis_range Integer vector of B-spline basis sizes, each at least 4.
#' @param topvarper Proportion in `(0, 1]` of genes used for tuning.
#' @param ncores Positive integer giving the number of parallel workers. Use
#'   `ncores = 1` for sequential execution.
#'
#' @return A list containing the selected `best_L`, `best_nbasis`, and
#'   `best_r_pen`, the minimum `GCV`, the vector `GCV_results`, and the evaluated
#'   `parameter_grid`.
#'
#' @examples
#' \donttest{
#' data(scFPCDE_simdata)
#' cells <- seq(1, 1000, length.out = 100)
#' cells <- unique(round(cells))
#' tuning <- scFPCDE_tune_fpca(
#'   scFPCDE_simdata$yt[cells, 1:40],
#'   scFPCDE_simdata$tt[cells],
#'   L = 2,
#'   r_pen_range = c(1e-3, 1e-2),
#'   nbasis_range = c(15, 20),
#'   ncores = 1
#' )
#' tuning[c("best_nbasis", "best_r_pen", "GCV")]
#' }
#' @export
scFPCDE_tune_fpca <- function(yt, tt, L = 3, r_pen_range = exp(-5:1),
                              nbasis_range = 1:5 * 20, topvarper = 0.1,
                              ncores = 2) {
  scFPCDE_validate_expression(yt, tt)
  L <- scFPCDE_validate_integer(L, "L")
  ncores <- scFPCDE_validate_integer(ncores, "ncores")
  scFPCDE_validate_probability(topvarper, "topvarper")
  if (!is.numeric(r_pen_range) || !length(r_pen_range) || anyNA(r_pen_range) ||
      any(!is.finite(r_pen_range)) || any(r_pen_range < 0)) {
    stop("`r_pen_range` must contain finite, non-negative numbers.", call. = FALSE)
  }
  if (!is.numeric(nbasis_range) || !length(nbasis_range) ||
      anyNA(nbasis_range) || any(!is.finite(nbasis_range)) ||
      any(nbasis_range < 4) || any(nbasis_range != as.integer(nbasis_range))) {
    stop("`nbasis_range` must contain integers greater than or equal to 4.",
         call. = FALSE)
  }
  if (any(L > nbasis_range)) {
    stop("`L` cannot exceed a value in `nbasis_range`.", call. = FALSE)
  }

  p <- ncol(yt)
  fpca0 <- scFPCDE_fit_fpca(
    yt, tt, L = L, r_pen = r_pen_range[1], nbasis = nbasis_range[1],
    fpc_varmax = FALSE
  )
  signal_strength <- rowSums(fpca0$scores^2)
  top_n <- ceiling(topvarper * p)
  top_idx <- utils::head(order(signal_strength, decreasing = TRUE), top_n)

  param_grid <- expand.grid(
    nbasis = as.integer(nbasis_range),
    r_pen = r_pen_range,
    KEEP.OUT.ATTRS = FALSE
  )

  calc_gcv <- function(params) {
    nbasis <- as.numeric(params["nbasis"])
    r_pen <- as.numeric(params["r_pen"])
    fpca_res <- scFPCDE_fit_fpca(
      yt = yt,
      tt = tt,
      L = L,
      r_pen = r_pen,
      nbasis = nbasis,
      topvarsub = top_idx,
      fpc_varmax = FALSE
    )
    mean(fpca_res$fda_splines$gcv)
  }

  if (ncores == 1L) {
    GCV_results <- vapply(
      seq_len(nrow(param_grid)),
      function(i) calc_gcv(param_grid[i, ]),
      numeric(1)
    )
  } else {
    cl <- scFPCDE_make_cluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    GCV_results <- parallel::parSapply(
      cl,
      seq_len(nrow(param_grid)),
      function(i) calc_gcv(param_grid[i, ])
    )
  }

  if (all(is.na(GCV_results))) {
    stop("All GCV evaluations returned missing values.", call. = FALSE)
  }
  best_idx <- which.min(GCV_results)
  best_params <- param_grid[best_idx, , drop = FALSE]

  list(
    best_L = L,
    best_nbasis = best_params$nbasis,
    best_r_pen = best_params$r_pen,
    GCV = min(GCV_results, na.rm = TRUE),
    GCV_results = GCV_results,
    parameter_grid = param_grid
  )
}
