#' Extract and Transform Gene-Level FPC Scores
#'
#' Extracts gene-level functional principal component (FPC) scores from an
#' scFPCDE fit or analysis result, optionally selects components, and applies a
#' visualization-oriented transformation.
#'
#' @param x A numeric score matrix, the result from [scFPCDE_fit_fpca()], or the
#'   result from [scFPCDE_run()]. Score matrices must have genes in rows and FPCs
#'   in columns.
#' @param components Optional FPCs to retain. Supply positive integer positions
#'   or component names. The default, `NULL`, retains every component.
#' @param transform Character string specifying the output scale: `"none"`
#'   returns raw scores, `"signed_log10"` applies a sign-preserving logarithmic
#'   compression, and `"standardize"` returns column-wise z-scores.
#' @param log_scale Positive number controlling the transition between the
#'   approximately linear and logarithmic regions of the signed-log transform.
#'   The default of 1 reproduces the transformation used in the manuscript
#'   analysis.
#'
#' @return A numeric matrix with genes in rows and selected FPCs in columns.
#'   The `"standardize"` result has `"scaled:center"` and `"scaled:scale"`
#'   attributes. Every result has a `"score_transform"` attribute, and a
#'   signed-log result also has a `"log_scale"` attribute.
#'
#' @details
#' The signed-log transformation is
#'
#' \deqn{\mathrm{sign}(x)\log_{10}(1 + |x|/s),}
#'
#' where \eqn{s} is `log_scale`. It preserves zero and the sign of each score,
#' and therefore preserves score-space quadrants while compressing extreme
#' coordinates. Standardization subtracts each component mean and divides by
#' its sample standard deviation. A constant component is centered and divided
#' by 1, producing zeros instead of missing values.
#'
#' These transformations are intended for visualization and exploratory
#' comparison. Differential-expression statistics and null boundaries should
#' be calculated from the raw scores returned by `transform = "none"`.
#' Component names produced by the underlying FPCA fit (`PC1`, `PC2`, and so
#' on) are normalized to the package terminology (`FPC1`, `FPC2`, and so on)
#' in the returned matrix; the fitted object itself is not modified.
#'
#' @examples
#' data(scFPCDE_simdata)
#' cells <- unique(round(seq(1, 1000, length.out = 80)))
#' fit <- scFPCDE_fit_fpca(
#'   scFPCDE_simdata$yt[cells, 1:20],
#'   scFPCDE_simdata$tt[cells],
#'   L = 2,
#'   nbasis = 12
#' )
#'
#' raw_scores <- scFPCDE_fpc_scores(fit)
#' compressed_scores <- scFPCDE_fpc_scores(
#'   fit,
#'   transform = "signed_log10"
#' )
#' standardized_scores <- scFPCDE_fpc_scores(
#'   fit,
#'   components = 1:2,
#'   transform = "standardize"
#' )
#' colMeans(standardized_scores)
#'
#' @export
scFPCDE_fpc_scores <- function(x, components = NULL,
                               transform = c("none", "signed_log10", "standardize"),
                               log_scale = 1) {
  transform <- match.arg(transform)
  scores <- scFPCDE_extract_scores(x)

  if (is.null(colnames(scores)) || anyNA(colnames(scores)) ||
      any(!nzchar(colnames(scores)))) {
    colnames(scores) <- paste0("FPC", seq_len(ncol(scores)))
  } else if (identical(colnames(scores), paste0("PC", seq_len(ncol(scores))))) {
    colnames(scores) <- paste0("FPC", seq_len(ncol(scores)))
  }

  if (!is.null(components)) {
    components <- scFPCDE_validate_components(components, scores)
    scores <- scores[, components, drop = FALSE]
  }

  if (transform == "signed_log10") {
    valid_log_scale <- length(log_scale) == 1L && is.numeric(log_scale) &&
      !is.na(log_scale) && is.finite(log_scale) && log_scale > 0
    if (!valid_log_scale) {
      stop("`log_scale` must be a single positive number.", call. = FALSE)
    }

    scores <- sign(scores) * log10(1 + abs(scores) / log_scale)
    attr(scores, "log_scale") <- log_scale
  } else if (transform == "standardize") {
    score_center <- colMeans(scores)
    score_scale <- apply(scores, 2L, stats::sd)
    score_scale[!is.finite(score_scale) | score_scale == 0] <- 1

    scores <- sweep(scores, 2L, score_center, "-")
    scores <- sweep(scores, 2L, score_scale, "/")
    attr(scores, "scaled:center") <- score_center
    attr(scores, "scaled:scale") <- score_scale
  }

  attr(scores, "score_transform") <- transform
  scores
}

scFPCDE_extract_scores <- function(x) {
  scores <- NULL

  if (is.matrix(x)) {
    scores <- x
  } else if (is.list(x) && is.matrix(x$scores)) {
    scores <- x$scores
  } else if (is.list(x) && is.list(x$fpca_result) &&
             is.matrix(x$fpca_result$scores)) {
    scores <- x$fpca_result$scores
  }

  if (is.null(scores)) {
    stop(
      "`x` must be a numeric score matrix or a result from `scFPCDE_fit_fpca()` or `scFPCDE_run()`.",
      call. = FALSE
    )
  }
  if (!is.numeric(scores) || !nrow(scores) || !ncol(scores)) {
    stop("The FPC score matrix must be non-empty and numeric.", call. = FALSE)
  }
  if (anyNA(scores) || any(!is.finite(scores))) {
    stop("The FPC score matrix must contain only finite, non-missing values.",
         call. = FALSE)
  }

  scores
}

scFPCDE_validate_components <- function(components, scores) {
  if (!length(components) || anyNA(components) || anyDuplicated(components)) {
    stop("`components` must select at least one unique FPC.", call. = FALSE)
  }

  if (is.numeric(components)) {
    valid <- all(is.finite(components)) && all(components %% 1 == 0) &&
      all(components >= 1) && all(components <= ncol(scores))
    if (!valid) {
      stop("Numeric `components` must be valid positive FPC positions.",
           call. = FALSE)
    }
    return(as.integer(components))
  }

  if (is.character(components)) {
    missing_components <- setdiff(components, colnames(scores))
    if (length(missing_components)) {
      stop(
        "Unknown FPC component(s): ",
        paste(missing_components, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    return(components)
  }

  stop("`components` must contain integer positions or component names.",
       call. = FALSE)
}
