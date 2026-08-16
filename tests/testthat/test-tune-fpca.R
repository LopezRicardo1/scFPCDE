test_that("FPCA tuning reports the best evaluated parameter pair", {
  dat <- small_scFPCDE_data(gene_index = 1:16)
  tuning <- scFPCDE_tune_fpca(
    dat$yt,
    dat$tt,
    L = 2,
    r_pen_range = c(1e-3, 1e-2),
    nbasis_range = c(10, 12),
    topvarper = 0.5,
    ncores = 1
  )

  expect_named(
    tuning,
    c(
      "best_L", "best_nbasis", "best_r_pen", "GCV", "GCV_results",
      "parameter_grid"
    )
  )
  expect_equal(nrow(tuning$parameter_grid), 4)
  expect_length(tuning$GCV_results, 4)
  expect_equal(tuning$GCV, min(tuning$GCV_results, na.rm = TRUE))
})
