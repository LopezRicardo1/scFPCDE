#' Compute Permutation P-Values
#'
#' Internal helper for computing permutation p-values using either a direct
#' comparison or the equivalent reverse-rank calculation.
#'
#' @param stats0 Numeric vector of observed statistics.
#' @param stats.perm Numeric vector containing the pooled permutation
#'   statistics.
#' @param method Character string selecting `"naive"` or `"revrank"`.
#'
#' @return A numeric vector of permutation p-values, one per observed
#'   statistic.
#' @keywords internal
permp <- function(stats0, stats.perm, method = c("naive", "revrank")) {
  method <- match.arg(method)
  if (!is.numeric(stats0) || !length(stats0) || anyNA(stats0)) {
    stop("`stats0` must be a non-empty numeric vector without missing values.",
         call. = FALSE)
  }
  if (!is.numeric(stats.perm) || !length(stats.perm) || anyNA(stats.perm)) {
    stop("`stats.perm` must be a non-empty numeric vector without missing values.",
         call. = FALSE)
  }

  m <- length(stats0)
  m.perm <- length(stats.perm)

  if (method == "naive") {
    pvals <- vapply(
      seq_len(m),
      function(i) sum(stats0[i] <= stats.perm) / m.perm,
      numeric(1)
    )
  } else {
    o <- order(stats0, decreasing = TRUE)
    stats.perm <- sort(stats.perm, decreasing = TRUE)
    j <- 1L
    rv <- numeric(m)
    for (i in seq_len(m)) {
      s0 <- stats0[o[i]]
      while (j <= m.perm && s0 <= stats.perm[j]) {
        j <- j + 1L
      }
      rv[o[i]] <- j - 1
    }
    pvals <- rv / m.perm
  }

  pvals
}
