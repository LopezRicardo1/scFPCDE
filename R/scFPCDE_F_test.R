#' FPCA F-statistic Permutation Test
#'
#' Performs an F-statistic based permutation test comparing null and FPCA models.
#'
#' @param yt A matrix of smoothed gene expression values (timepoints × genes).
#' @param tt Time vector.
#' @param fpca_result A list with FPCA results: `$scores`, `$eigenfunctions`, and `$sigma2`.
#' @param n_perm Number of permutations.
#' @param ncores Number of parallel cores to use for permutation. Default is 2.
#' @return A data frame with gene IDs, observed F statistic, p-values, and q-values.
#' @export
#' @importFrom stats p.adjust
scFPCDE_F_test <- function(yt, tt, fpca_result, n_perm = 1000, ncores = 2) {
  scores <- t(fpca_result$scores)
  Phi <- fpca_result$eigenfunctions
  sigma2 <- fpca_result$sigma2

  RSS0 <- colSums(yt^2)
  yt_hat <- Phi %*% scores
  RSS1 <- colSums((yt - yt_hat)^2)
  F_obs <- (RSS0 - RSS1) / (RSS1 + sigma2)

  cl <- scFPCDE_make_cluster(ncores = ncores)
  permuted_F_stats <- parallel::parSapply(cl, 1:n_perm, function(i) {
    yt_perm <- yt[sample(nrow(yt)), ]
    scores_perm <- solve(crossprod(Phi), crossprod(Phi, yt_perm))
    yt_perm_hat <- Phi %*% scores_perm
    RSS1_perm <- colSums((yt_perm - yt_perm_hat)^2)
    (RSS0 - RSS1_perm) / (RSS1_perm + sigma2)
  })
  parallel::stopCluster(cl)

  p_val <- permp(F_obs, as.vector(permuted_F_stats))
  q_val <- p.adjust(p_val, "BH")

  data.frame(
    ID = colnames(yt),
    F_obs = F_obs,
    p_value = p_val,
    q_value = q_val
  )
}
