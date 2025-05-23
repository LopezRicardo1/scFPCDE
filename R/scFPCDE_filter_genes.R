#' Gene Filtering by Expression Proportion
#'
#' Filters genes based on the proportion of non-zero expression values across cells.
#'
#' @param y A gene expression matrix (cells × genes)
#' @param qz Quantile threshold for filtering (default = 0.25)
#' @return A character vector of filtered gene names
#' @export
scFPCDE_filter_genes <- function(y, qz = 0.25) {
  z <- ifelse(y == 0, 0, 1)
  sumz <- colSums(z)
  threshold <- quantile(sumz, qz)
  filtered_genes <- colnames(y[, sumz > threshold])
  return(filtered_genes)
}
