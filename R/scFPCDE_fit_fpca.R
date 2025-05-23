#' Fit Functional Principal Component Analysis (FPCA)
#'
#' Performs FPCA on smoothed gene expression curves using B-spline basis.
#'
#' @param yt Matrix of gene expression (rows = timepoints, columns = genes)
#' @param tt Time vector
#' @param L Number of principal components (harmonics)
#' @param r_pen Smoothing penalty
#' @param nbasis Number of basis functions
#' @param topvarsub Optional indices of genes to subset for FPCA basis estimation
#'
#' @return A list with scores, eigenfunctions, fitted values, and FPCA internals
#' @export
scFPCDE_fit_fpca <- function(yt, tt, L = 2, r_pen = 1e-3, nbasis = 50, topvarsub = NULL) {
  original_yt <- yt
  if (!is.null(topvarsub)) yt <- yt[, topvarsub]

  yt_centered <- sweep(yt, 2, colMeans(yt), "-")
  basis <- fda::create.bspline.basis(range(tt), nbasis = nbasis)
  par <- fda::fdPar(basis, 2, lambda = r_pen)
  ss <- fda::smooth.basis(tt, yt_centered, par)
  fpca <- fda::pca.fd(ss$fd, nharm = L, centerfns = FALSE)

  Phi <- fda::eval.fd(tt, fpca$harmonics)
  original_yt_centered <- sweep(original_yt, 2, colMeans(original_yt), "-")
  scores <- solve(crossprod(Phi), crossprod(Phi, original_yt_centered))
  xt_hat <- Phi %*% scores
  sigma2 <- mean((original_yt_centered - xt_hat)^2)

  return(list(
    scores = t(scores),
    eigenfunctions = Phi,
    eigenvals = apply(scores, 1, var),
    PEV = fpca$varprop,
    sigma2 = sigma2,
    xt_hat = xt_hat,
    fda_splines = ss,
    fda_fpca = fpca
  ))
}
