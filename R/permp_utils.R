#' Reverse-Rank Permutation P-values
#'
#' Helper function to compute permutation p-values using either naive or reverse-rank method.
#'
#' @param stats0 Vector of observed statistics
#' @param stats.perm Vector of permutation statistics
#' @param method Method to compute p-values: "naive" or "revrank"
#' @return Vector of p-values
#' @keywords internal
permp <- function(stats0, stats.perm, method = c("naive", "revrank")) {
  method <- match.arg(method)
  m <- length(stats0); m.perm <- length(stats.perm)

  if (method == "naive") {
    pvals <- sapply(1:m, function(i) sum(stats0[i] <= stats.perm)) / m.perm
  } else {
    o <- order(stats0, decreasing = TRUE)
    stats.perm <- sort(stats.perm, decreasing = TRUE)
    j <- 1; rv <- rep(0, m)
    for (i in 1:m) {
      s0 <- stats0[o[i]]
      while (s0 <= stats.perm[j] && j <= m.perm) j <- j + 1
      rv[o[i]] <- j - 1
    }
    pvals <- rv / m.perm
  }

  return(pvals)
}
