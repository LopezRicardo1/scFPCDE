#' Simulated scFPCDE Dataset
#'
#' A simulated dataset with 1000 cells and 500 genes, including pseudotime and cluster labels.
#' 100 genes are differentially expressed using functional basis combinations.
#' 400 are null genes.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{yt}{A numeric matrix of gene expression (1000 cells × 500 genes)}
#'   \item{tt}{A numeric vector of pseudotime values (length 1000)}
#'   \item{clusters}{A character vector of cluster assignments ("1", "2", "3")}
#' }
#' @usage data(scFPCDE_simdata)
#' @source Generated with `scFPCDE_simulate_data()` (now internal)
"scFPCDE_simdata"
