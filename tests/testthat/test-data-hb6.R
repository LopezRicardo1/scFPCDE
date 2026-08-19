test_that("HB6 trajectories preserve the paper-analysis dimensions", {
  data(scFPCDE_hb6, package = "scFPCDE", envir = environment())

  expect_named(scFPCDE_hb6, c("traj1", "traj2"))
  expect_equal(dim(scFPCDE_hb6$traj1$yt), c(2928, 2995))
  expect_equal(dim(scFPCDE_hb6$traj2$yt), c(1714, 1999))
  expect_identical(
    levels(scFPCDE_hb6$traj1$clusters),
    c("Trans", "Naive", "C-mem1")
  )
  expect_identical(
    levels(scFPCDE_hb6$traj2$clusters),
    c("M-mem1", "DN2", "DN3")
  )
})

test_that("HB6 expression, pseudotime, and labels are exactly aligned", {
  data(scFPCDE_hb6, package = "scFPCDE", envir = environment())

  for (trajectory in scFPCDE_hb6) {
    expect_identical(rownames(trajectory$yt), trajectory$cell_id)
    expect_identical(colnames(trajectory$yt), trajectory$gene_id)
    expect_identical(dimnames(trajectory$counts), dimnames(trajectory$yt))
    expect_identical(names(trajectory$tt), trajectory$cell_id)
    expect_identical(names(trajectory$clusters), trajectory$cell_id)
    expect_equal(anyDuplicated(trajectory$gene_id), 0L)
    expect_true(all(is.finite(trajectory$yt)))
    expect_true(all(is.finite(trajectory$counts)))
    expect_true(all(trajectory$counts >= 0))
    expect_true(all(trajectory$counts == floor(trajectory$counts)))
    expect_true(all(is.finite(trajectory$tt)))
    expect_false(anyNA(trajectory$clusters))
    expect_false(isTRUE(all.equal(
      unname(colMeans(trajectory$yt)),
      rep(0, ncol(trajectory$yt))
    )))
    expect_false(trajectory$metadata$expression_centered)
    expect_identical(trajectory$metadata$expression_assay, "logcounts")
    expect_identical(trajectory$metadata$count_assay, "counts")
  }
})

test_that("HB6 provenance records the public study and paper settings", {
  data(scFPCDE_hb6, package = "scFPCDE", envir = environment())

  for (trajectory in scFPCDE_hb6) {
    expect_identical(trajectory$metadata$study$donor, "HB6")
    expect_identical(trajectory$metadata$study$accession, "E-MTAB-9544")
    expect_identical(
      trajectory$metadata$study$doi,
      "10.3389/fimmu.2021.602539"
    )
    expect_equal(trajectory$metadata$paper_parameters$L, 3)
    expect_equal(trajectory$metadata$paper_parameters$nbasis, 30)
    expect_equal(trajectory$metadata$paper_parameters$n_perm, 500)
  }
})
