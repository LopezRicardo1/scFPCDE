#' Cross-Platform Cluster Creation Helper
#'
#' Creates a parallel cluster using `makeForkCluster()` on Unix-like systems and
#' `makeCluster()` on Windows for compatibility.
#'
#' @param ncores Number of cores to use (default = detectCores() - 1)
#' @return A cluster object from the parallel package
#' @keywords internal
scFPCDE_make_cluster <- function(ncores = ncores <- min(2, max(1, parallel::detectCores() - 1))) {
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(ncores)
  } else {
    cl <- parallel::makeForkCluster(ncores)
  }
  return(cl)
}
