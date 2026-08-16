#' Run the FPCA F-Statistic Permutation Test
#'
#' Compares residual sums of squares from the null and fitted FPCA models and
#' evaluates the resulting statistic against a pooled permutation distribution.
#'
#' @param yt Numeric centered expression matrix with observations in rows and
#'   genes in columns.
#' @param tt Numeric pseudotime vector with one value per row of `yt`.
#' @param fpca_result An FPCA result returned by [scFPCDE_fit_fpca()].
#' @param n_perm Positive integer giving the number of permutations.
#' @param ncores Positive integer giving the number of parallel workers. Use
#'   `ncores = 1` for sequential execution.
#'
#' @return A data frame with one row per gene and columns `ID`, `F_obs`,
#'   `p_value`, and `q_value`. The q-values use the Benjamini-Hochberg method.
#'
#' @details
#' The function uses R's current random-number state. For deterministic results,
#' call [set.seed()] and use `ncores = 1`.
#'
#' @examples
#' \donttest{
#' data(scFPCDE_simdata)
#' cells <- seq(1, 1000, length.out = 100)
#' cells <- unique(round(cells))
#' yt <- scale(scFPCDE_simdata$yt[cells, 1:40], center = TRUE, scale = FALSE)
#' tt <- scFPCDE_simdata$tt[cells]
#' fit <- scFPCDE_fit_fpca(yt, tt, nbasis = 20)
#' set.seed(1)
#' test <- scFPCDE_F_test(yt, tt, fit, n_perm = 10, ncores = 1)
#' head(test)
#' }
#' @export
scFPCDE_F_test <- function(yt, tt, fpca_result, n_perm = 1000, ncores = 2) {
  scFPCDE_validate_expression(yt, tt)
  scFPCDE_validate_fpca_result(fpca_result, yt, require_sigma2 = TRUE)
  n_perm <- scFPCDE_validate_integer(n_perm, "n_perm")
  ncores <- scFPCDE_validate_integer(ncores, "ncores")

  scores <- t(fpca_result$scores)
  Phi <- fpca_result$eigenfunctions
  sigma2 <- fpca_result$sigma2

  RSS0 <- colSums(yt^2)
  yt_hat <- Phi %*% scores
  RSS1 <- colSums((yt - yt_hat)^2)
  F_obs <- (RSS0 - RSS1) / (RSS1 + sigma2)

  calculate_permutation <- function(i) {
    yt_perm <- yt[sample.int(nrow(yt)), , drop = FALSE]
    scores_perm <- solve(crossprod(Phi), crossprod(Phi, yt_perm))
    yt_perm_hat <- Phi %*% scores_perm
    RSS1_perm <- colSums((yt_perm - yt_perm_hat)^2)
    (RSS0 - RSS1_perm) / (RSS1_perm + sigma2)
  }

  if (ncores == 1L) {
    permuted_F_stats <- vapply(
      seq_len(n_perm), calculate_permutation, numeric(ncol(yt))
    )
  } else {
    cl <- scFPCDE_make_cluster(ncores = ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    permuted_F_stats <- parallel::parSapply(
      cl, seq_len(n_perm), calculate_permutation
    )
  }

  p_val <- permp(F_obs, as.vector(permuted_F_stats))
  q_val <- stats::p.adjust(p_val, "BH")

  data.frame(
    ID = colnames(yt),
    F_obs = F_obs,
    p_value = p_val,
    q_value = q_val
  )
}
