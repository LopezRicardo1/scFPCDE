#' Plot Observed and Fitted Gene-Expression Curves
#'
#' Displays FPCA-fitted trajectories over observed expression values. Points are
#' colored by cell cluster, and genes can be faceted or overlaid in one panel.
#'
#' @param tt Numeric pseudotime vector with one value per row of `yt`.
#' @param yt Numeric observed expression matrix with observations in rows and
#'   genes in columns. When plotting a result from [scFPCDE_run()], apply the
#'   same centering and scaling used in that call.
#' @param yt_fit Numeric matrix of fitted expression values with the same
#'   dimensions and gene names as `yt`.
#' @param cell_cluster Vector of cluster assignments with one value per row of
#'   `yt`.
#' @param subset Optional vector of gene names or indices to plot. Use `NULL` to
#'   plot every gene.
#' @param facet_genes Logical; if `TRUE`, draw genes in separate facets.
#' @param nrow Optional number of facet rows.
#' @param ncol Optional number of facet columns.
#' @param legend_dot_size Positive numeric size for points in the cluster legend.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @examples
#' \donttest{
#' data(scFPCDE_simdata)
#' cells <- seq(1, 1000, length.out = 100)
#' cells <- unique(round(cells))
#' yt <- scale(scFPCDE_simdata$yt[cells, 1:12], center = TRUE, scale = FALSE)
#' fit <- scFPCDE_fit_fpca(yt, scFPCDE_simdata$tt[cells], nbasis = 20)
#' scFPCDE_gene_curves(
#'   tt = scFPCDE_simdata$tt[cells],
#'   yt = yt,
#'   yt_fit = fit$xt_hat,
#'   cell_cluster = scFPCDE_simdata$clusters[cells]
#' )
#' }
#' @export
#' @importFrom magrittr %>%
scFPCDE_gene_curves <- function(tt, yt, yt_fit, cell_cluster,
                                     subset = 1:12,
                                     facet_genes = TRUE,
                                     nrow = 4, ncol = 3,
                                     legend_dot_size = 3) {

  scFPCDE_validate_expression(yt, tt)
  if (!is.matrix(yt_fit) || !is.numeric(yt_fit) ||
      !identical(dim(yt_fit), dim(yt)) || anyNA(yt_fit) ||
      any(!is.finite(yt_fit))) {
    stop("`yt_fit` must be a finite numeric matrix with the same dimensions as `yt`.",
         call. = FALSE)
  }
  if (!identical(colnames(yt_fit), colnames(yt))) {
    stop("`yt_fit` and `yt` must have identical gene names.", call. = FALSE)
  }
  if (length(cell_cluster) != nrow(yt) || anyNA(cell_cluster)) {
    stop("`cell_cluster` must contain one non-missing value per row of `yt`.",
         call. = FALSE)
  }
  if (length(facet_genes) != 1L || is.na(facet_genes) || !is.logical(facet_genes)) {
    stop("`facet_genes` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (length(legend_dot_size) != 1L || !is.numeric(legend_dot_size) ||
      is.na(legend_dot_size) || !is.finite(legend_dot_size) ||
      legend_dot_size <= 0) {
    stop("`legend_dot_size` must be a single positive number.", call. = FALSE)
  }

  if (!is.null(subset)) {
    if (!length(subset) || anyNA(subset)) {
      stop("`subset` must select at least one gene.", call. = FALSE)
    }
    selected <- tryCatch(
      list(
        yt = yt[, subset, drop = FALSE],
        yt_fit = yt_fit[, subset, drop = FALSE]
      ),
      error = function(e) {
        stop("`subset` contains invalid gene indices or names.", call. = FALSE)
      }
    )
    yt <- selected$yt
    yt_fit <- selected$yt_fit
    if (!ncol(yt)) {
      stop("`subset` must select at least one gene.", call. = FALSE)
    }
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
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = legend_dot_size))) +
    ggplot2::labs(y = "Expression", title = "Gene Expression Functions by Cluster") +
    ggplot2::facet_wrap(~Gene, nrow = nrow, ncol = ncol, scales = "free_y")

  p_color <- ggplot2::ggplot(combined_data, ggplot2::aes(Pseudotime, Fitted_Expression, color = Gene)) +
    ggplot2::geom_point(ggplot2::aes(y = Expression, colour = Cluster), size = 0.01, alpha = 0.25) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme(aspect.ratio = .6) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = legend_dot_size))) +
    ggplot2::labs(y = "Expression", title = "Gene Expression Functions")

  if (facet_genes) p_facet else p_color
}
