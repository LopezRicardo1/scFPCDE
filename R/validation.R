scFPCDE_validate_expression <- function(yt, tt = NULL, argument = "yt") {
  if (!is.matrix(yt) || !is.numeric(yt)) {
    stop("`", argument, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(yt) < 2L || ncol(yt) < 1L) {
    stop("`", argument, "` must contain at least two rows and one column.", call. = FALSE)
  }
  if (anyNA(yt) || any(!is.finite(yt))) {
    stop("`", argument, "` must contain only finite, non-missing values.", call. = FALSE)
  }
  if (is.null(colnames(yt))) {
    stop("`", argument, "` must have gene names as column names.", call. = FALSE)
  }

  if (!is.null(tt)) {
    if (!is.numeric(tt) || length(tt) != nrow(yt)) {
      stop("`tt` must be a numeric vector with one value per row of `", argument, "`.", call. = FALSE)
    }
    if (anyNA(tt) || any(!is.finite(tt))) {
      stop("`tt` must contain only finite, non-missing values.", call. = FALSE)
    }
    if (length(unique(tt)) < 2L) {
      stop("`tt` must contain at least two distinct pseudotime values.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

scFPCDE_validate_integer <- function(x, argument, minimum = 1L) {
  valid <- length(x) == 1L && is.numeric(x) && !is.na(x) && is.finite(x)
  if (valid) {
    valid <- x >= minimum && x <= .Machine$integer.max && x %% 1 == 0
  }
  if (!valid) {
    stop("`", argument, "` must be a single integer greater than or equal to ",
         minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

scFPCDE_validate_probability <- function(x, argument, include_zero = FALSE) {
  valid <- length(x) == 1L && is.numeric(x) && !is.na(x) && is.finite(x)
  if (valid) {
    lower_ok <- if (include_zero) x >= 0 else x > 0
    valid <- lower_ok && x <= 1
  }
  if (!valid) {
    interval <- if (include_zero) "[0, 1]" else "(0, 1]"
    stop("`", argument, "` must be a single number in ", interval, ".", call. = FALSE)
  }
  invisible(TRUE)
}

scFPCDE_validate_fpca_arguments <- function(L, r_pen, nbasis, fpc_varmax) {
  L <- scFPCDE_validate_integer(L, "L")
  nbasis <- scFPCDE_validate_integer(nbasis, "nbasis", minimum = 4L)
  if (L > nbasis) {
    stop("`L` cannot exceed `nbasis`.", call. = FALSE)
  }
  if (length(r_pen) != 1L || is.na(r_pen) || !is.numeric(r_pen) ||
      !is.finite(r_pen) || r_pen < 0) {
    stop("`r_pen` must be a single non-negative number.", call. = FALSE)
  }
  if (length(fpc_varmax) != 1L || is.na(fpc_varmax) || !is.logical(fpc_varmax)) {
    stop("`fpc_varmax` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  invisible(TRUE)
}

scFPCDE_validate_fpca_result <- function(fpca_result, yt, require_sigma2 = FALSE) {
  required <- c("scores", "eigenfunctions")
  if (require_sigma2) {
    required <- c(required, "sigma2")
  }
  if (!is.list(fpca_result) || !all(required %in% names(fpca_result))) {
    stop("`fpca_result` must contain ", paste(required, collapse = ", "), ".",
         call. = FALSE)
  }

  scores <- fpca_result$scores
  eigenfunctions <- fpca_result$eigenfunctions
  if (!is.matrix(scores) || !is.numeric(scores) || nrow(scores) != ncol(yt) ||
      anyNA(scores) || any(!is.finite(scores))) {
    stop("`fpca_result$scores` must have one row per gene in `yt`.", call. = FALSE)
  }
  if (!is.matrix(eigenfunctions) || !is.numeric(eigenfunctions) ||
      nrow(eigenfunctions) != nrow(yt) || ncol(eigenfunctions) != ncol(scores)) {
    stop("`fpca_result$eigenfunctions` has incompatible dimensions.", call. = FALSE)
  }
  if (anyNA(eigenfunctions) || any(!is.finite(eigenfunctions))) {
    stop("`fpca_result$eigenfunctions` must contain only finite values.",
         call. = FALSE)
  }
  if (require_sigma2 &&
      (length(fpca_result$sigma2) != 1L || !is.numeric(fpca_result$sigma2) ||
       is.na(fpca_result$sigma2) || !is.finite(fpca_result$sigma2) ||
       fpca_result$sigma2 < 0)) {
    stop("`fpca_result$sigma2` must be a single non-negative number.", call. = FALSE)
  }

  invisible(TRUE)
}
