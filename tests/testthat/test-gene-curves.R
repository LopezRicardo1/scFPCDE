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
})

test_that("gene curve plotting requires one cluster label per observation", {
  dat <- small_scFPCDE_data(gene_index = 1:4)
  yt <- scale(dat$yt, center = TRUE, scale = FALSE)

  expect_error(
    scFPCDE_gene_curves(dat$tt, yt, yt, dat$clusters[-1], subset = 1:2),
    "one non-missing value per row"
  )
})
