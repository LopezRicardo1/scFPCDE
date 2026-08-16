test_that("FPCA fit preserves the documented dimensions", {
  dat <- small_scFPCDE_data()
  fit <- scFPCDE_fit_fpca(
    dat$yt,
    dat$tt,
    L = 2,
    nbasis = 12,
    fpc_varmax = FALSE
  )

  expect_named(
    fit,
    c(
      "scores", "eigenfunctions", "eigenvals", "PEV", "sigma2",
      "xt_hat", "fda_splines", "fda_fpca", "fpc_varmax"
    )
  )
  expect_equal(dim(fit$scores), c(ncol(dat$yt), 2))
  expect_equal(dim(fit$eigenfunctions), c(nrow(dat$yt), 2))
  expect_equal(dim(fit$xt_hat), dim(dat$yt))
  expect_identical(colnames(fit$xt_hat), colnames(dat$yt))
})

test_that("the wrapper returns stable, fully named result components", {
  dat <- small_scFPCDE_data(gene_index = c(1:8, 101:108))
  set.seed(42)
  result <- scFPCDE_run(
    dat$yt,
    dat$tt,
    L = 2,
    nbasis = 12,
    n_perm = 3,
    topvarper = 0.5,
    ncores = 1,
    use_FPC_F = TRUE,
    fpc_varmax = FALSE
  )

  expect_named(result, c("fpca_result", "D_test_result", "F_test_result"))
  expect_s3_class(result$D_test_result, "data.frame")
  expect_s3_class(result$F_test_result, "data.frame")
  expect_named(result$D_test_result, c("ID", "D_obs", "p_value", "q_value"))
  expect_named(result$F_test_result, c("ID", "F_obs", "p_value", "q_value"))
  expect_identical(result$D_test_result$ID, colnames(dat$yt))
  expect_true(all(result$D_test_result$p_value >= 0))
  expect_true(all(result$D_test_result$p_value <= 1))
})

test_that("sequential permutations are reproducible after set.seed", {
  dat <- small_scFPCDE_data()
  yt <- scale(dat$yt, center = TRUE, scale = FALSE)
  fit <- scFPCDE_fit_fpca(
    yt, dat$tt, L = 2, nbasis = 12, fpc_varmax = FALSE
  )

  set.seed(99)
  first <- scFPCDE_D_test(yt, dat$tt, fit, n_perm = 3, ncores = 1)
  set.seed(99)
  second <- scFPCDE_D_test(yt, dat$tt, fit, n_perm = 3, ncores = 1)

  expect_identical(first, second)
})
