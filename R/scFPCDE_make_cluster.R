#' Create a Cross-Platform Parallel Cluster
#'
#' Internal helper that creates a fork cluster on Unix-like systems and a PSOCK
#' cluster on Windows.
#'
#' @param ncores Positive integer giving the number of workers.
#'
#' @return A cluster object from [parallel::makeCluster()].
#' @keywords internal
scFPCDE_make_cluster <- function(ncores = 2) {
  ncores <- scFPCDE_validate_integer(ncores, "ncores")
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(ncores)
  } else {
    cl <- parallel::makeForkCluster(ncores)
  }

  cl
}
