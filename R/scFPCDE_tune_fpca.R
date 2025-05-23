#' Tune FPCA Parameters via Generalized Cross-Validation (GCV)
#'
#' Optimizes smoothing parameters for FPCA using GCV on top signal curves.
#'
#' @param yt Gene expression matrix (timepoints × genes)
#' @param tt Time vector
#' @param L Number of harmonics to retain
#' @param r_pen_range Vector of smoothing penalties
#' @param nbasis_range Vector of number of basis functions
#' @param topvarper Proportion of top signal genes used for tuning
#'
#' @return A list with optimal parameters and full GCV results
#' @export
scFPCDE_tune_fpca <- function(yt, tt, L = 3, r_pen_range = exp(-5:1),
                              nbasis_range = 1:5 * 20, topvarper = 0.1) {
  p <- ncol(yt)
  fpca0 <- scFPCDE_fit_fpca(yt, tt, L, r_pen = r_pen_range[1], nbasis = nbasis_range[1])
  signal_strength <- rowSums(fpca0$scores^2)
  top_n <- ceiling(topvarper * p)
  top_idx <- order(signal_strength, decreasing = TRUE)[1:top_n]

  param_grid <- expand.grid(nbasis = nbasis_range, r_pen = r_pen_range)

  calc_gcv <- function(params) {
    nbasis <- as.numeric(params["nbasis"])
    r_pen <- as.numeric(params["r_pen"])
    fpca_res <- scFPCDE_fit_fpca(yt, tt, L, r_pen, nbasis, topvarsub = top_idx)
    sum(fpca_res$fda_splines$gcv)
  }

  cl <- scFPCDE_make_cluster()
  GCV_results <- parallel::parSapply(cl, 1:nrow(param_grid), function(i) calc_gcv(param_grid[i, ]))
  parallel::stopCluster(cl)

  best_idx <- which(GCV_results == min(GCV_results, na.rm = TRUE))
  best_params <- param_grid[best_idx, ]

  return(list(
    best_L = L,
    best_nbasis = as.numeric(best_params["nbasis"]),
    best_r_pen = as.numeric(best_params["r_pen"]),
    GCV = min(GCV_results, na.rm = TRUE),
    GCV_results = GCV_results
  ))
}
