#' Cross-Platform Cluster Creation Helper
#'
#' Creates a parallel cluster using `makeForkCluster()` on Unix-like systems and
#' `makeCluster()` on Windows for compatibility.
#'
#' @param ncores Number of cores to use. Default is 2.
#' @return A cluster object from the parallel package.
#' @keywords internal
scFPCDE_make_cluster <- function(ncores = 2) {
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(ncores)
  } else {
    cl <- parallel::makeForkCluster(ncores)
  }

  return(cl)
}
