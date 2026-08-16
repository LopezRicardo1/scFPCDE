#' Fit Functional Principal Component Analysis
#'
#' Smooths centered gene-expression trajectories with a B-spline basis, fits
#' functional principal components, and projects every input gene onto the
#' estimated eigenfunctions.
#'
#' @param yt Numeric expression matrix with cells or time points in rows and
#'   genes in columns. Gene names are required as column names.
#' @param tt Numeric pseudotime vector with one value per row of `yt`.
#' @param L Positive integer giving the number of harmonics to retain.
#' @param r_pen Non-negative smoothing penalty passed to [fda::fdPar()].
#' @param nbasis Integer number of B-spline basis functions; must be at least 4.
#' @param topvarsub Optional gene indices or names used to estimate the FPCA
#'   eigenfunctions. All genes are still projected onto the fitted basis.
#' @param fpc_varmax Logical; if `TRUE`, apply VARIMAX rotation to the fitted
#'   harmonics.
#'
#' @return A list with components:
#' \describe{
#'   \item{scores}{A genes-by-`L` matrix of FPCA scores.}
#'   \item{eigenfunctions}{An observations-by-`L` matrix of evaluated
#'   eigenfunctions.}
#'   \item{eigenvals}{The variance of each score dimension across genes.}
#'   \item{PEV}{The proportion of variation explained by each harmonic.}
#'   \item{sigma2}{Mean squared residual variation across all genes.}
#'   \item{xt_hat}{The fitted centered expression matrix.}
#'   \item{fda_splines}{The object returned by [fda::smooth.basis()].}
#'   \item{fda_fpca}{The fitted object returned by [fda::pca.fd()].}
#'   \item{fpc_varmax}{Whether VARIMAX rotation was applied.}
#' }
#'
#' @examples
#' \donttest{
#' data(scFPCDE_simdata)
#' cells <- seq(1, 1000, length.out = 100)
#' cells <- unique(round(cells))
#' fit <- scFPCDE_fit_fpca(
#'   scFPCDE_simdata$yt[cells, 1:40],
#'   scFPCDE_simdata$tt[cells],
#'   L = 2,
#'   nbasis = 20
#' )
#' dim(fit$scores)
#' }
#' @export
scFPCDE_fit_fpca <- function(yt, tt, L = 2, r_pen = 1e-3, nbasis = 50, topvarsub = NULL, fpc_varmax = TRUE) {
  scFPCDE_validate_expression(yt, tt)
  scFPCDE_validate_fpca_arguments(L, r_pen, nbasis, fpc_varmax)

  original_yt <- yt
  if (!is.null(topvarsub)) {
    if (!length(topvarsub) || anyNA(topvarsub)) {
      stop("`topvarsub` must select at least one gene.", call. = FALSE)
    }
    yt <- tryCatch(
      yt[, topvarsub, drop = FALSE],
      error = function(e) {
        stop("`topvarsub` contains invalid gene indices or names.", call. = FALSE)
      }
    )
    if (!ncol(yt)) {
      stop("`topvarsub` must select at least one gene.", call. = FALSE)
    }
  }

  yt_centered <- sweep(yt, 2, colMeans(yt), "-")
  basis <- fda::create.bspline.basis(range(tt), nbasis = nbasis)
  par <- fda::fdPar(basis, 2, lambda = r_pen)
  ss <- fda::smooth.basis(tt, yt_centered, par)
  fpca <- fda::pca.fd(ss$fd, nharm = L, centerfns = FALSE)
  if (fpc_varmax) {
    fpca <- fda::varmx.pca.fd(fpca)
  }
  Phi <- fda::eval.fd(tt, fpca$harmonics)
  original_yt_centered <- sweep(original_yt, 2, colMeans(original_yt), "-")
  scores <- solve(crossprod(Phi), crossprod(Phi, original_yt_centered))
  xt_hat <- Phi %*% scores
  sigma2 <- mean((original_yt_centered - xt_hat)^2)

  list(
    scores = t(scores),
    eigenfunctions = Phi,
    eigenvals = apply(scores, 1, stats::var),
    PEV = fpca$varprop,
    sigma2 = sigma2,
    xt_hat = xt_hat,
    fda_splines = ss,
    fda_fpca = fpca,
    fpc_varmax = fpc_varmax
  )
}
