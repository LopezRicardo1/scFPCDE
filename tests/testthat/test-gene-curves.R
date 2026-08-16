test_that("gene curve plotting returns a ggplot object", {
  dat <- small_scFPCDE_data(gene_index = 1:6)
  yt <- scale(dat$yt, center = TRUE, scale = FALSE)
  fit <- scFPCDE_fit_fpca(
    yt, dat$tt, L = 2, nbasis = 12, fpc_varmax = FALSE
  )

  plot <- scFPCDE_gene_curves(
    dat$tt,
    yt,
    fit$xt_hat,
    dat$clusters,
    subset = 1:4,
    nrow = 2,
    ncol = 2
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$x, "Pseudotime")
  expect_equal(plot$labels$y, "Centered expression")
  expect_equal(plot$labels$colour, "Cell cluster")
  expect_equal(nlevels(plot$data$Gene), 4L)
})

test_that("gene curve plotting requires one cluster label per observation", {
  dat <- small_scFPCDE_data(gene_index = 1:4)
  yt <- scale(dat$yt, center = TRUE, scale = FALSE)

  expect_error(
    scFPCDE_gene_curves(dat$tt, yt, yt, dat$clusters[-1], subset = 1:2),
    "one non-missing value per row"
  )
})

test_that("gene curve plotting supports cohesive visual controls", {
  dat <- small_scFPCDE_data(gene_index = 1:6)
  yt <- scale(dat$yt, center = TRUE, scale = FALSE)
  fit <- scFPCDE_fit_fpca(
    yt, dat$tt, L = 2, nbasis = 12, fpc_varmax = FALSE
  )
  cluster_levels <- unique(as.character(dat$clusters))
  colors <- stats::setNames(
    c("#0072B2", "#009E73", "#D55E00")[seq_along(cluster_levels)],
    cluster_levels
  )

  plot <- scFPCDE_gene_curves(
    dat$tt,
    yt,
    fit$xt_hat,
    dat$clusters,
    subset = 1:3,
    facet_genes = FALSE,
    cluster_colors = colors,
    point_size = 0.5,
    point_alpha = 0.6,
    curve_linewidth = 1.1,
    title = "Selected trajectories"
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$title, "Selected trajectories")
  expect_equal(plot$labels$linetype, "Gene")
  expect_equal(
    unname(plot$scales$get_scales("colour")$palette(length(colors))),
    unname(colors)
  )
})

test_that("gene curve plotting validates visual controls", {
  dat <- small_scFPCDE_data(gene_index = 1:4)
  yt <- scale(dat$yt, center = TRUE, scale = FALSE)

  expect_error(
    scFPCDE_gene_curves(
      dat$tt, yt, yt, dat$clusters, subset = 1:2,
      point_alpha = 1.5
    ),
    "point_alpha"
  )
  expect_error(
    scFPCDE_gene_curves(
      dat$tt, yt, yt, dat$clusters, subset = 1:2,
      cluster_colors = "black"
    ),
    "one color per cluster"
  )
})
