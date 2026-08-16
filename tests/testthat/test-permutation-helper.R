test_that("reverse-rank and direct permutation p-values agree", {
  observed <- c(3, 2, 1)
  permuted <- c(4, 2.5, 2, 0)

  direct <- scFPCDE:::permp(observed, permuted, method = "naive")
  reverse_rank <- scFPCDE:::permp(observed, permuted, method = "revrank")

  expect_equal(reverse_rank, direct)
})

test_that("reverse-rank handles a statistic below every permutation", {
  expect_equal(
    scFPCDE:::permp(0, c(1, 2, 3), method = "revrank"),
    1
  )
})
