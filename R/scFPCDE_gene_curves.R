#' Plot Observed and Fitted Gene-Expression Curves
#'
#' Displays FPCA-fitted trajectories over observed expression values. Observed
#' cells are colored by cluster and fitted trajectories are drawn as cohesive
#' curves. Genes can be faceted or overlaid in one panel.
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
#' @param point_size Positive numeric size for observed-expression points.
#' @param point_alpha Numeric transparency for observed-expression points in
#'   `[0, 1]`.
#' @param curve_linewidth Positive numeric width for fitted trajectories.
#' @param curve_color Color used for fitted trajectories.
#' @param cluster_colors Optional vector of colors for cell clusters. A named
#'   vector is matched to cluster labels; an unnamed vector is used in cluster
#'   level order.
#' @param scales Facet scale behavior passed to [ggplot2::facet_wrap()]. One of
#'   `"free_y"`, `"fixed"`, `"free"`, or `"free_x"`.
#' @param title Optional plot title. Use `NULL` for a context-specific default.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @details
#' All rows of `tt`, `yt`, `yt_fit`, and `cell_cluster` must describe the same
#' observations in the same order. The function orders these inputs together by
#' pseudotime before plotting. Duplicate gene names are disambiguated in facet
#' labels with [make.unique()].
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
#'   cell_cluster = scFPCDE_simdata$clusters[cells],
#'   point_size = 0.5,
#'   point_alpha = 0.6
#' )
#' }
#' @export
scFPCDE_gene_curves <- function(tt, yt, yt_fit, cell_cluster,
                                subset = 1:12,
                                facet_genes = TRUE,
                                nrow = 4, ncol = 3,
                                legend_dot_size = 3,
                                point_size = 0.45,
                                point_alpha = 0.55,
                                curve_linewidth = 0.9,
                                curve_color = "black",
                                cluster_colors = NULL,
                                scales = "free_y",
                                title = NULL) {

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
  for (argument in c("point_size", "curve_linewidth")) {
    value <- get(argument)
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value <= 0) {
      stop("`", argument, "` must be a single positive number.", call. = FALSE)
    }
  }
  scFPCDE_validate_probability(point_alpha, "point_alpha", include_zero = TRUE)
  if (length(curve_color) != 1L || !is.character(curve_color) ||
      is.na(curve_color) || !nzchar(curve_color)) {
    stop("`curve_color` must be a single color string.", call. = FALSE)
  }
  scales <- match.arg(scales, c("free_y", "fixed", "free", "free_x"))
  if (!is.null(title) &&
      (length(title) != 1L || !is.character(title) || is.na(title))) {
    stop("`title` must be `NULL` or a single character string.", call. = FALSE)
  }
  for (argument in c("nrow", "ncol")) {
    value <- get(argument)
    if (!is.null(value)) {
      value <- scFPCDE_validate_integer(value, argument)
      if (argument == "nrow") nrow <- value else ncol <- value
    }
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

  ord <- order(tt)
  tt <- tt[ord]
  yt <- yt[ord, , drop = FALSE]
  yt_fit <- yt_fit[ord, , drop = FALSE]
  cell_cluster <- cell_cluster[ord]

  cluster <- if (is.factor(cell_cluster)) {
    droplevels(cell_cluster)
  } else {
    factor(cell_cluster, levels = unique(as.character(cell_cluster)))
  }
  cluster_levels <- levels(cluster)

  if (!is.null(cluster_colors)) {
    if ((!is.character(cluster_colors) && !is.numeric(cluster_colors)) ||
        !length(cluster_colors) || anyNA(cluster_colors)) {
      stop("`cluster_colors` must be a non-missing vector of colors.",
           call. = FALSE)
    }
    if (is.null(names(cluster_colors))) {
      if (length(cluster_colors) < length(cluster_levels)) {
        stop("`cluster_colors` must provide one color per cluster.",
             call. = FALSE)
      }
      cluster_colors <- stats::setNames(
        cluster_colors[seq_along(cluster_levels)], cluster_levels
      )
    } else {
      missing_colors <- setdiff(cluster_levels, names(cluster_colors))
      if (length(missing_colors)) {
        stop(
          "`cluster_colors` is missing colors for: ",
          paste(missing_colors, collapse = ", "), ".",
          call. = FALSE
        )
      }
      cluster_colors <- cluster_colors[cluster_levels]
    }
  }

  gene_labels <- make.unique(colnames(yt))
  combined_data <- data.frame(
    Pseudotime = rep(tt, times = ncol(yt)),
    Cluster = rep(cluster, times = ncol(yt)),
    Gene = factor(
      rep(gene_labels, each = nrow(yt)),
      levels = gene_labels
    ),
    Expression = as.vector(yt),
    Fitted_Expression = as.vector(yt_fit),
    check.names = FALSE
  )

  if (is.null(title)) {
    title <- if (facet_genes) {
      "Observed and FPCA-fitted gene trajectories"
    } else {
      "FPCA-fitted gene trajectories"
    }
  }

  plot <- ggplot2::ggplot(combined_data, ggplot2::aes(x = Pseudotime)) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", linewidth = 0.35,
      color = "grey65"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = Expression, colour = Cluster),
      size = point_size, alpha = point_alpha
    ) +
    ggplot2::labs(
      x = "Pseudotime",
      y = "Centered expression",
      colour = "Cell cluster",
      title = title
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = legend_dot_size, alpha = 1)
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey92", linewidth = 0.25
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey96", color = "grey75"
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(cluster_colors)) {
    plot <- plot + ggplot2::scale_colour_manual(values = cluster_colors)
  }

  if (facet_genes) {
    plot <- plot +
      ggplot2::geom_line(
        ggplot2::aes(y = Fitted_Expression, group = Gene),
        color = curve_color,
        linewidth = curve_linewidth,
        lineend = "round"
      ) +
      ggplot2::facet_wrap(
        ~Gene, nrow = nrow, ncol = ncol, scales = scales
      )
  } else {
    plot <- plot +
      ggplot2::geom_line(
        ggplot2::aes(
          y = Fitted_Expression,
          group = Gene,
          linetype = Gene
        ),
        color = curve_color,
        linewidth = curve_linewidth,
        lineend = "round"
      ) +
      ggplot2::labs(linetype = "Gene")
  }

  plot
}
