#' Run Full scFPCDE Wrapper Pipeline
#'
#' Wrapper function to run the full scFPCDE workflow: smoothing with FPCA,
#' basis refinement via top-variance filtering, and permutation tests.
#'
#' @param yt Matrix of gene expression (rows = timepoints, columns = genes).
#' @param tt Time vector.
#' @param L Number of harmonics (latent dimensions).
#' @param r_pen Smoothing penalty.
#' @param nbasis Number of basis functions.
#' @param n_perm Number of permutations.
#' @param topvarper Proportion of top variable genes for eigenfunction estimation.
#' @param center Logical; center each gene across time (default = TRUE).
#' @param scale Logical; scale each gene to unit variance across time (default = FALSE).
#' @param ncores Number of parallel cores to use for permutation tests. Default is 2.
#' @param use_FPC_F Logical; whether to run the FPCA-based F-test. Default is FALSE.
#'
#' @return A list with FPCA result, D-test, and optionally F-test results.
#' @export
scFPCDE_run <- function(yt, tt, L = 2, r_pen = 1e-3, nbasis = 50,
                        n_perm = 1000, topvarper = 0.1,
                        center = TRUE, scale = FALSE,
                        ncores = 2, use_FPC_F = FALSE) {

  # Order by pseudotime to ensure consistent trajectory alignment
  ord <- order(tt)
  tt <- tt[ord]
  yt <- yt[ord, , drop = FALSE]

  # Preprocess gene expression: center and/or scale each gene
  yt <- scale(yt, center = center, scale = scale)

  # Initial FPCA using all genes
  fpca_result_full <- scFPCDE_fit_fpca(
    yt, tt, L = L, r_pen = r_pen, nbasis = nbasis
  )

  # Select top variable genes
  D_stat <- rowSums(fpca_result_full$scores^2)
  n_genes <- ncol(yt)
  top_n <- ceiling(topvarper * n_genes)
  top_idx <- order(D_stat, decreasing = TRUE)[1:top_n]

  # Refit FPCA on top variable genes
  fpca_result_top <- scFPCDE_fit_fpca(
    yt, tt, L = L, r_pen = r_pen, nbasis = nbasis,
    topvarsub = top_idx
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

  return(list(
    fpca_result = fpca_result_top,
    D_test_result = D_test_result,
    F_test_result = F_test_result
  ))
}
