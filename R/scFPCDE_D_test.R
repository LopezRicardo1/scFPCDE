#' FPCA Distance-Based Permutation Test (D-statistic)
#'
#' Performs a permutation-based D-statistic test using functional PCA scores.
#'
#' @param yt A matrix of smoothed gene expression values (timepoints × genes)
#' @param tt Time vector
#' @param fpca_result A list with FPCA results: `$scores` and `$eigenfunctions`
#' @param n_perm Number of permutations
#' @return A data frame with gene IDs, observed D statistic, p-values, and q-values
#' @export
#' @importFrom stats p.adjust
scFPCDE_D_test <- function(yt, tt, fpca_result, n_perm = 1000) {
  scores <- fpca_result$scores
  Phi <- fpca_result$eigenfunctions
  D_obs <- rowSums(scores^2)

  cl <- scFPCDE_make_cluster()
  permuted_D_stats <- parallel::parSapply(cl, 1:n_perm, function(i) {
    yt_perm <- yt[sample(nrow(yt)), ]
    scores_perm <- t(solve(crossprod(Phi), crossprod(Phi, yt_perm)))
    rowSums(scores_perm^2)
  })
  parallel::stopCluster(cl)

  p_val <- permp(D_obs, as.vector(permuted_D_stats))
  q_val <- p.adjust(p_val, "BH")

  data.frame(
    ID = colnames(yt),
    D_obs = D_obs,
    p_value = p_val,
    q_value = q_val
  )
}
