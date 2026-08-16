test_that("gene filtering uses the requested detection-count quantile", {
  y <- cbind(
    gene_a = c(1, 0, 0, 0),
    gene_b = c(1, 1, 0, 0),
    gene_c = c(1, 1, 1, 1)
  )

  expect_equal(scFPCDE_filter_genes(y, qz = 0.25), c("gene_b", "gene_c"))
  expect_identical(scFPCDE_filter_genes(y, qz = 1), character())
})

test_that("gene filtering validates its public inputs", {
  y <- matrix(1:6, nrow = 2)

  expect_error(scFPCDE_filter_genes(y), "gene names")
  colnames(y) <- paste0("gene_", 1:3)
  expect_error(scFPCDE_filter_genes(y, qz = -0.1), "\\[0, 1\\]")
})

test_that("duplicate gene names remain backward compatible", {
  y <- cbind(
    gene_a = c(1, 0, 0, 0),
    gene_a = c(1, 1, 0, 0),
    gene_b = c(1, 1, 1, 1)
  )

  expect_equal(scFPCDE_filter_genes(y, qz = 0.25), c("gene_a", "gene_b"))
})
