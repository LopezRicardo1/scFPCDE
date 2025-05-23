#' Plot Smoothed and Observed Gene Expression Curves
#'
#' Visualizes fitted gene expression trajectories from FPCA alongside observed values, colored by cluster.
#'
#' @param tt Vector of pseudotime values (length = # of timepoints)
#' @param yt Matrix of original expression data (rows = timepoints, columns = genes)
#' @param yt_fit Matrix of fitted expression values (same dimensions as yt)
#' @param cell_cluster Vector of cluster assignments (length = ncol(yt))
#' @param subset Optional vector of gene names or indices to subset
#' @param facet_genes If TRUE, plots genes in separate facets (default = TRUE)
#' @param nrow Number of facet rows (if faceted)
#' @param ncol Number of facet columns (if faceted)
#' @param legend_dot_size Size of legend dots for cluster colors
#'
#' @return A ggplot object
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats var quantile
scFPCDE_gene_curves <- function(tt, yt, yt_fit, cell_cluster,
                                     subset = NULL,
                                     facet_genes = TRUE,
                                     nrow = 7, ncol = 10,
                                     legend_dot_size = 3) {

  if (!is.null(subset)) {
    yt <- yt[, subset, drop = FALSE]
    yt_fit <- yt_fit[, subset, drop = FALSE]
  }

  y_long <- data.frame(Pseudotime = tt, yt, check.names = FALSE) %>%
    tidyr::pivot_longer(cols = -Pseudotime, names_to = "Gene", values_to = "Expression", names_repair = "minimal")

  y_fit_long <- data.frame(Pseudotime = tt, Cluster = cell_cluster, yt_fit, check.names = FALSE) %>%
    tidyr::pivot_longer(cols = -(1:2), names_to = "Gene", values_to = "Fitted_Expression", names_repair = "minimal")

  combined_data <- dplyr::mutate(y_fit_long, Expression = y_long$Expression)

  if (!is.null(subset)) {
    combined_data$Gene <- factor(combined_data$Gene, levels = colnames(yt))
  }

  p_facet <- ggplot2::ggplot(combined_data, ggplot2::aes(Pseudotime, Fitted_Expression)) +
    ggplot2::geom_point(ggplot2::aes(y = Expression, colour = Cluster), size = 0.01, alpha = 0.25) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = legend_dot_size))) +
    ggplot2::labs(y = "Expression", title = "Gene Expression Functions by Cluster") +
    ggplot2::facet_wrap(~Gene, nrow = nrow, ncol = ncol, scales = "free_y")

  p_color <- ggplot2::ggplot(combined_data, ggplot2::aes(Pseudotime, Fitted_Expression, color = Gene)) +
    ggplot2::geom_point(ggplot2::aes(y = Expression, colour = Cluster), size = 0.01, alpha = 0.25) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme(aspect.ratio = .6) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = legend_dot_size))) +
    ggplot2::labs(y = "Expression", title = "Gene Expression Functions")

  return(if (facet_genes) p_facet else p_color)
}
