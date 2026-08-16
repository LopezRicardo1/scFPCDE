small_scFPCDE_data <- function(n_cells = 80L, gene_index = 1:12) {
  data(scFPCDE_simdata, package = "scFPCDE", envir = environment())
  cell_index <- unique(round(seq(1, 1000, length.out = n_cells)))

  list(
    yt = scFPCDE_simdata$yt[cell_index, gene_index, drop = FALSE],
    tt = scFPCDE_simdata$tt[cell_index],
    clusters = scFPCDE_simdata$clusters[cell_index]
  )
}
