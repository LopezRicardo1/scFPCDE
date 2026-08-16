test_that("FPC scores are extracted from supported package results", {
  scores <- matrix(
    c(-9, -1, 0, 1, 9, 2, 2, 2, 2, 2),
    ncol = 2,
    dimnames = list(paste0("gene", 1:5), c("FPC1", "FPC2"))
  )

  raw_matrix <- scFPCDE_fpc_scores(scores)
  raw_fit <- scFPCDE_fpc_scores(list(scores = scores))
  raw_run <- scFPCDE_fpc_scores(list(fpca_result = list(scores = scores)))

  expect_equal(raw_matrix, scores, ignore_attr = TRUE)
  expect_equal(raw_fit, scores, ignore_attr = TRUE)
  expect_equal(raw_run, scores, ignore_attr = TRUE)
  expect_identical(dimnames(raw_matrix), dimnames(scores))
  expect_identical(attr(raw_matrix, "score_transform"), "none")
})

test_that("FPC scores integrate with fitted package objects", {
  dat <- small_scFPCDE_data(n_cells = 50, gene_index = 1:8)
  fit <- scFPCDE_fit_fpca(
    dat$yt,
    dat$tt,
    L = 2,
    nbasis = 8,
    fpc_varmax = FALSE
  )

  scores <- scFPCDE_fpc_scores(
    fit,
    components = c("FPC2", "FPC1"),
    transform = "signed_log10"
  )

  expect_identical(dim(scores), c(8L, 2L))
  expect_identical(rownames(scores), colnames(dat$yt))
  expect_identical(colnames(scores), c("FPC2", "FPC1"))
  expect_true(all(is.finite(scores)))
})

test_that("signed-log compression matches the manuscript transformation", {
  scores <- matrix(
    c(-99, -9, 0, 9, 99, -3, -1, 0, 1, 3),
    ncol = 2,
    dimnames = list(paste0("gene", 1:5), c("FPC1", "FPC2"))
  )
  original <- scores

  transformed <- scFPCDE_fpc_scores(
    scores,
    transform = "signed_log10"
  )

  expected <- sign(scores) * log10(1 + abs(scores))
  expect_equal(unname(transformed), unname(expected), ignore_attr = TRUE)
  expect_equal(sign(transformed), sign(scores), ignore_attr = TRUE)
  expect_identical(attr(transformed, "log_scale"), 1)
  expect_identical(attr(transformed, "score_transform"), "signed_log10")
  expect_identical(scores, original)

  scaled <- scFPCDE_fpc_scores(
    scores,
    transform = "signed_log10",
    log_scale = 10
  )
  expect_equal(
    unname(scaled),
    unname(sign(scores) * log10(1 + abs(scores) / 10)),
    ignore_attr = TRUE
  )
})

test_that("standardization produces stable column-wise z-scores", {
  scores <- cbind(
    FPC1 = c(-2, -1, 0, 1, 2),
    FPC2 = rep(4, 5)
  )
  rownames(scores) <- paste0("gene", 1:5)

  standardized <- scFPCDE_fpc_scores(
    scores,
    transform = "standardize"
  )

  expect_equal(mean(standardized[, "FPC1"]), 0)
  expect_equal(stats::sd(standardized[, "FPC1"]), 1)
  expect_equal(unname(standardized[, "FPC2"]), rep(0, 5), ignore_attr = TRUE)
  expect_equal(attr(standardized, "scaled:center"), c(FPC1 = 0, FPC2 = 4))
  expect_equal(
    attr(standardized, "scaled:scale"),
    c(FPC1 = stats::sd(scores[, "FPC1"]), FPC2 = 1)
  )
  expect_identical(attr(standardized, "score_transform"), "standardize")
  expect_identical(dimnames(standardized), dimnames(scores))
})

test_that("FPC components can be selected by position or name", {
  scores <- matrix(
    seq_len(18),
    nrow = 6,
    dimnames = list(paste0("gene", 1:6), paste0("FPC", 1:3))
  )

  selected_position <- scFPCDE_fpc_scores(scores, components = c(3, 1))
  selected_name <- scFPCDE_fpc_scores(
    scores,
    components = c("FPC2", "FPC1")
  )
  expect_equal(selected_position, scores[, c(3, 1), drop = FALSE], ignore_attr = TRUE)
  expect_identical(colnames(selected_position), c("FPC3", "FPC1"))
  expect_equal(
    selected_name,
    scores[, c("FPC2", "FPC1"), drop = FALSE],
    ignore_attr = TRUE
  )
  expect_identical(colnames(selected_name), c("FPC2", "FPC1"))

  unnamed_scores <- unname(scores)
  selected <- scFPCDE_fpc_scores(unnamed_scores, components = "FPC2")
  expect_identical(colnames(selected), "FPC2")
})

test_that("invalid score transformations fail clearly", {
  scores <- matrix(seq_len(12), nrow = 4)

  expect_error(scFPCDE_fpc_scores(data.frame(scores)), "must be a numeric score matrix")
  expect_error(scFPCDE_fpc_scores(matrix(character(), 0, 0)), "non-empty and numeric")
  expect_error(
    scFPCDE_fpc_scores(matrix(c(1, NA), nrow = 1)),
    "finite, non-missing"
  )
  expect_error(scFPCDE_fpc_scores(scores, components = integer()), "at least one")
  expect_error(scFPCDE_fpc_scores(scores, components = c(1, 1)), "unique FPC")
  expect_error(scFPCDE_fpc_scores(scores, components = 4), "valid positive")
  expect_error(scFPCDE_fpc_scores(scores, components = "FPC9"), "Unknown FPC")
  expect_error(
    scFPCDE_fpc_scores(scores, transform = "signed_log10", log_scale = 0),
    "single positive"
  )
})
