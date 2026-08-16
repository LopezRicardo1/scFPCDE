#' Filter Genes by Detection Frequency
#'
#' Retains genes whose number of non-zero observations is greater than a
#' user-specified quantile of the gene-level detection counts.
#'
#' @param y Numeric expression matrix with cells in rows and genes in columns.
#'   Column names are required and are returned for genes that pass the filter.
#' @param qz Numeric quantile in `[0, 1]` used as the detection-count threshold.
#'
#' @return A character vector containing the names of genes with detection
#'   counts strictly greater than the selected quantile.
#'
#' @examples
#' y <- matrix(
#'   c(0, 0, 1, 1, 1, 1, 0, 1, 2, 3, 1, 2),
#'   nrow = 4,
#'   dimnames = list(NULL, c("gene_a", "gene_b", "gene_c"))
#' )
#' scFPCDE_filter_genes(y, qz = 0.25)
#' @export
scFPCDE_filter_genes <- function(y, qz = 0.25) {
  scFPCDE_validate_expression(y, argument = "y")
  scFPCDE_validate_probability(qz, "qz", include_zero = TRUE)

  z <- ifelse(y == 0, 0, 1)
  sumz <- colSums(z)
  threshold <- stats::quantile(sumz, qz)

  colnames(y)[sumz > threshold]
}
