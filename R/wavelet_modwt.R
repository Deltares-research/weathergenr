# ==============================================================================
# Maximal Overlap Discrete Wavelet Transform (MODWT)
# ==============================================================================
#
# Shift-invariant wavelet decomposition for time series analysis.
# Provides additive multi-resolution decomposition at dyadic scales.
#
# This module wraps the waveslim package for robust MODWT/MRA computation.
#
# Main functions:
# - modwt_decompose(): Forward MODWT decomposition
# - modwt_reconstruct(): Inverse MODWT reconstruction
# - modwt_mra(): Multiresolution analysis (additive components)
# - analyze_wavelet_additive(): Hybrid CWT + MODWT analysis
#
# Terminology note:
# MODWT MRA provides ADDITIVE (not orthogonal) decomposition. The components
# sum exactly to the original series, but are NOT statistically independent.
# Cross-covariances between components are generally non-zero.
#
# Low-frequency safeguard:
# By default, the maximum represented period is capped to <= (1/3)*n (record length),
# to reduce spurious very-low-frequency structure in short records.
#
# References:
# Percival, D. B., & Walden, A. T. (2000). Wavelet Methods for Time Series Analysis.
# Cambridge University Press.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Wavelet helpers (internal)
# ------------------------------------------------------------------------------

#' Get wavelet filter length
#'
#' @description
#' Returns the filter length for a named wavelet filter used by MODWT/MRA routines.
#'
#' @param filter Character scalar. Filter name.
#'
#' @return Integer scalar. Filter length.
#'
#' @keywords internal
.get_filter_length <- function(filter) {
  lengths <- c(
    haar = 2L,
    d4   = 4L,
    d6   = 6L,
    d8   = 8L,
    la8  = 8L,
    la16 = 16L
  )

  if (!is.character(filter) || length(filter) != 1L || is.na(filter)) {
    stop("'filter' must be a non-missing character scalar.", call. = FALSE)
  }

  if (!filter %in% names(lengths)) {
    stop("Unknown filter: ", filter, call. = FALSE)
  }

  lengths[[filter]]
}

#' Validate boundary handling option
#'
#' @param boundary Character scalar. Boundary option.
#' @return Character scalar. Normalized boundary.
#' @keywords internal
.check_boundary <- function(boundary) {
  if (!is.character(boundary) || length(boundary) != 1L || is.na(boundary)) {
    stop("'boundary' must be a non-missing character scalar.", call. = FALSE)
  }
  if (!identical(boundary, "periodic")) {
    stop("Only 'periodic' boundary is supported.", call. = FALSE)
  }
  boundary
}

#' Validate a numeric series input
#'
#' @param x Vector. Input series.
#' @param name Character scalar. Argument name (for error messages).
#' @param min_n Integer scalar. Minimum required length.
#' @return Numeric vector. Coerced to numeric.
#' @keywords internal
.as_numeric_series <- function(x, name, min_n) {
  if (!is.numeric(x)) stop("'", name, "' must be a numeric vector.", call. = FALSE)
  if (anyNA(x)) stop("'", name, "' contains missing values.", call. = FALSE)

  x <- as.numeric(x)
  n <- length(x)

  if (!is.finite(min_n) || length(min_n) != 1L) {
    stop("Internal error: 'min_n' must be a scalar.", call. = FALSE)
  }
  if (n < as.integer(min_n)) {
    stop("'", name, "' must have at least ", as.integer(min_n), " observations.", call. = FALSE)
  }

  x
}

#' Validate a positive integer scalar argument
#'
#' @param x Object. Candidate scalar.
#' @param name Character scalar. Argument name.
#' @return Integer scalar.
#' @keywords internal
.as_pos_int_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1) {
    stop("'", name, "' must be a positive integer.", call. = FALSE)
  }
  as.integer(x)
}

#' Compute max allowed MODWT levels given a max-period fraction
#'
#' @description
#' Ensures the smooth threshold period 2^(J+1) does not exceed floor(n * max_period_frac).
#' This prevents representing variability at periods longer than a chosen fraction of the
#' record length.
#'
#' @param n Integer scalar. Series length.
#' @param max_period_frac Numeric scalar in (0, 1]. Maximum allowed period as a fraction of n.
#'
#' @return Integer scalar. Maximum allowed number of levels (J), at least 1.
#'
#' @keywords internal
.max_levels_by_period_frac <- function(n, max_period_frac) {
  if (!is.numeric(max_period_frac) || length(max_period_frac) != 1L ||
      !is.finite(max_period_frac) || max_period_frac <= 0 || max_period_frac > 1) {
    stop("'max_period_frac' must be in (0, 1].", call. = FALSE)
  }

  n <- .as_pos_int_scalar(n, "n")

  max_period_allowed <- floor(as.numeric(n) * max_period_frac)

  # Need at least 4 to allow J >= 1 because 2^(J+1) with J = 1 gives 4
  if (!is.finite(max_period_allowed) || max_period_allowed < 4L) return(1L)

  Jmax <- floor(log2(max_period_allowed)) - 1L
  max(1L, as.integer(Jmax))
}

#' Compute recommended maximum MODWT levels based on length and filter width
#'
#' @param n Integer scalar. Series length.
#' @param filter_length Integer scalar. Wavelet filter length.
#' @param max_period_frac Numeric scalar in (0, 1]. Optional cap based on record-length fraction.
#' @return Integer scalar. Recommended maximum levels (J), at least 1.
#' @keywords internal
.recommended_max_levels <- function(n, filter_length, max_period_frac = 1) {
  n <- .as_pos_int_scalar(n, "n")
  filter_length <- .as_pos_int_scalar(filter_length, "filter_length")

  max_levels_length <- floor(log2(n))
  max_levels_filter <- floor(log2((n - 1) / (filter_length - 1) + 1))
  max_levels <- max(1L, min(as.integer(max_levels_length), as.integer(max_levels_filter)))

  max_levels_period <- .max_levels_by_period_frac(n, max_period_frac)
  min(max_levels, max_levels_period)
}

#' Choose n_levels given requested and recommended maximum
#'
#' @param n_levels Integer scalar or NULL.
#' @param max_levels Integer scalar.
#' @param default_backoff Integer scalar. If NULL n_levels, use max(1, max_levels - default_backoff).
#' @return Integer scalar. Selected n_levels in [1, max_levels].
#' @keywords internal
.select_n_levels <- function(n_levels, max_levels, default_backoff = 1L) {
  max_levels <- .as_pos_int_scalar(max_levels, "max_levels")

  if (is.null(n_levels)) {
    out <- max(1L, max_levels - as.integer(default_backoff))
    return(as.integer(out))
  }

  n_levels <- .as_pos_int_scalar(n_levels, "n_levels")

  if (n_levels > max_levels) {
    warning(
      sprintf(
        "Requested n_levels=%d exceeds recommended maximum=%d. Using %d level(s).",
        n_levels, max_levels, max_levels
      ),
      call. = FALSE
    )
    return(as.integer(max_levels))
  }

  as.integer(n_levels)
}

# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------

#' Maximal Overlap Discrete Wavelet Transform (MODWT)
#'
#' @description
#' Performs the Maximal Overlap Discrete Wavelet Transform using a specified wavelet filter.
#' MODWT is shift-invariant and provides additive multi-resolution decomposition at dyadic
#' scales (powers of 2).
#'
#' This function wraps \code{waveslim::modwt()} with input validation and a standardized
#' return object.
#'
#' @param x Numeric vector. Input time series. Must be regularly spaced and contain no missing
#'   values.
#' @param filter Character. Wavelet filter name. One of \code{"la8"} (default), \code{"haar"},
#'   \code{"d4"}, \code{"d6"}, \code{"d8"}, \code{"la16"}.
#' @param n_levels Integer scalar or NULL. Number of decomposition levels (J). If NULL, a
#'   conservative default is used based on series length and filter width.
#' @param boundary Character. Boundary handling method. Only \code{"periodic"} is supported.
#'
#' @return A list with class \code{"modwt_result"}.
#' The object contains:
#' - \code{coefficients}: coefficients returned by \code{waveslim::modwt()} (W1...WJ and VJ)
#' - \code{filter}: filter name used
#' - \code{n_levels}: number of decomposition levels (J)
#' - \code{n}: series length
#' - \code{boundary}: boundary method
#' - \code{filter_length}: length of the wavelet filter
#'
#' @details
#' The recommended maximum number of levels is capped using a stability heuristic based on
#' series length and filter width. The cap is intended to reduce boundary artifacts and
#' unstable behavior on short records; it is not a mathematical requirement of MODWT.
#'
#' @importFrom waveslim modwt
#' @keywords internal
#' @export
modwt_decompose <- function(
    x,
    filter = c("la8", "haar", "d4", "d6", "d8", "la16"),
    n_levels = NULL,
    boundary = "periodic"
) {
  x <- .as_numeric_series(x, "x", min_n = 8L)
  boundary <- .check_boundary(boundary)
  filter <- match.arg(filter)

  n <- length(x)
  filter_length <- .get_filter_length(filter)
  max_levels <- .recommended_max_levels(n = n, filter_length = filter_length, max_period_frac = 1)

  # default: back off by 1 for a conservative choice
  n_levels <- .select_n_levels(n_levels, max_levels, default_backoff = 1L)

  wt <- waveslim::modwt(x, wf = filter, n.levels = n_levels, boundary = boundary)

  structure(
    list(
      coefficients = wt,
      filter = filter,
      n_levels = n_levels,
      n = n,
      boundary = boundary,
      filter_length = filter_length
    ),
    class = "modwt_result"
  )
}

#' Inverse MODWT
#'
#' @description
#' Reconstructs a time series from MODWT coefficients.
#'
#' @param wt Either a \code{"modwt_result"} object returned by \code{\link{modwt_decompose}},
#'   or a raw \code{waveslim::modwt()} output list containing W* and VJ (or V<level>) entries.
#'
#' @return Numeric vector. Reconstructed time series.
#'
#' @importFrom waveslim imodwt
#' @keywords internal
#' @export
modwt_reconstruct <- function(wt) {
  if (inherits(wt, "modwt_result")) {
    coeffs <- wt$coefficients

  } else if (is.list(wt)) {
    nm <- names(wt)
    has_W <- any(grepl("^W[0-9]+$", nm))
    has_V <- any(grepl("^V[0-9]+$", nm)) || ("VJ" %in% nm)

    if (!has_W || !has_V) {
      stop(
        "'wt' must be a 'modwt_result' object or valid waveslim modwt output ",
        "(requires W* and VJ or V<level> coefficients).",
        call. = FALSE
      )
    }

    coeffs <- wt

  } else {
    stop("'wt' must be a 'modwt_result' object or waveslim modwt output.", call. = FALSE)
  }

  as.numeric(waveslim::imodwt(coeffs))
}

#' MODWT multiresolution analysis (MRA)
#'
#' @description
#' Extracts additive time-domain components from a MODWT decomposition using multiresolution
#' analysis. Each component captures variance within a dyadic scale band. The components
#' sum (approximately, under periodic boundary handling) to the original series.
#'
#' @param x Numeric vector. Input time series. Must contain no missing values.
#' @param filter Character. Wavelet filter name. One of \code{"la8"} (default), \code{"haar"},
#'   \code{"d4"}, \code{"d6"}, \code{"d8"}, \code{"la16"}.
#' @param n_levels Integer scalar or NULL. Number of decomposition levels (J). If NULL, a
#'   conservative default is selected subject to stability and max-period caps.
#' @param boundary Character. Boundary handling method. Only \code{"periodic"} is supported.
#' @param include_smooth Logical. If TRUE, includes the smooth (approximation) component at
#'   the coarsest level (SJ).
#' @param max_period_frac Numeric scalar in (0, 1]. Maximum represented period as a fraction
#'   of record length. Default is \code{1/3}, meaning no structure beyond approximately n/3.
#'
#' @return A list with additive components and summary diagnostics. The list contains:
#' - \code{components}: numeric matrix (n x n_components) with columns D1..DJ and optionally SJ
#' - \code{periods}: numeric vector of representative periods per component
#' - \code{variance}: component variances
#' - \code{variance_fraction}: component variance divided by total variance of x
#' - \code{total_variance}: variance of x
#' - \code{sum_component_variance}: sum of component variances
#' - \code{cross_covariance_sum}: sum of cross-covariances (pairwise, without doubles)
#' - \code{variance_identity_error}: total_var - (sum_var + 2 * cross_cov_sum)
#' - \code{reconstruction_error}: max absolute error between x and reconstructed series
#' - \code{n_levels}, \code{filter}, \code{filter_length}, \code{boundary}, \code{include_smooth},
#'   \code{max_period_frac}, and \code{is_additive}
#'
#' @importFrom waveslim mra
#' @export
modwt_mra <- function(
    x,
    filter = c("la8", "haar", "d4", "d6", "d8", "la16"),
    n_levels = NULL,
    boundary = "periodic",
    include_smooth = TRUE,
    max_period_frac = 1/3
) {
  x <- .as_numeric_series(x, "x", min_n = 8L)
  boundary <- .check_boundary(boundary)
  filter <- match.arg(filter)

  n <- length(x)
  filter_length <- .get_filter_length(filter)

  max_levels <- .recommended_max_levels(
    n = n,
    filter_length = filter_length,
    max_period_frac = max_period_frac
  )

  n_levels <- .select_n_levels(n_levels, max_levels, default_backoff = 1L)

  mra_result <- waveslim::mra(
    x,
    wf = filter,
    J = n_levels,
    method = "modwt",
    boundary = boundary
  )

  # Components: D1..DJ and optionally SJ
  n_comp <- if (isTRUE(include_smooth)) n_levels + 1L else n_levels
  components <- matrix(0, nrow = n, ncol = n_comp)

  detail_names <- paste0("D", seq_len(n_levels))
  smooth_name <- paste0("S", n_levels)

  for (j in seq_len(n_levels)) {
    components[, j] <- mra_result[[detail_names[j]]]
  }

  if (isTRUE(include_smooth)) {
    components[, n_comp] <- mra_result[[smooth_name]]
  }

  col_names <- c(detail_names, if (isTRUE(include_smooth)) smooth_name else NULL)
  colnames(components) <- col_names

  # Representative periods
  # - Detail j roughly corresponds to band [2^j, 2^(j+1)] -> report geometric mean 2^j * sqrt(2)
  # - Smooth threshold period is 2^(J+1)
  detail_periods <- 2^seq_len(n_levels) * sqrt(2)
  smooth_threshold <- 2^(n_levels + 1)
  periods <- c(detail_periods, if (isTRUE(include_smooth)) smooth_threshold else NULL)

  # Variance accounting
  comp_var <- apply(components, 2, stats::var)
  total_var <- stats::var(x)
  sum_comp_var <- sum(comp_var)

  cov_mat <- stats::cov(components)
  cross_cov <- cov_mat
  diag(cross_cov) <- 0
  cross_cov_sum <- sum(cross_cov) / 2

  variance_identity_error <- total_var - (sum_comp_var + 2 * cross_cov_sum)

  var_frac <- comp_var / total_var

  reconstructed <- rowSums(components)
  recon_error <- max(abs(x - reconstructed))

  list(
    components = components,
    periods = periods,
    variance = comp_var,
    variance_fraction = var_frac,
    total_variance = total_var,
    sum_component_variance = sum_comp_var,
    cross_covariance_sum = cross_cov_sum,
    variance_identity_error = variance_identity_error,
    n_levels = n_levels,
    filter = filter,
    filter_length = filter_length,
    boundary = boundary,
    include_smooth = isTRUE(include_smooth),
    max_period_frac = max_period_frac,
    reconstruction_error = recon_error,
    is_additive = TRUE
  )
}

#' Additive wavelet analysis with CWT diagnostics
#'
#' @description
#' Combines CWT analysis (visualization and significance testing) with MODWT MRA (additive
#' components). CWT is used for diagnostics only; MODWT MRA provides the additive component
#' extraction used downstream.
#'
#' This function expects \code{analyze_wavelet_spectrum()} to be available elsewhere in the
#' package.
#'
#' @param series Numeric vector. Input series (regularly spaced, no missing values).
#' @param signif Numeric scalar in (0, 1). Significance level for CWT.
#' @param noise Character. Background noise model for CWT significance testing. Typically
#'   \code{"white"} or \code{"red"}.
#' @param min_period Numeric scalar. Minimum Fourier period for CWT.
#' @param detrend Logical. If TRUE, removes a linear trend before CWT.
#' @param filter Character. MODWT filter name. One of \code{"la8"} (default), \code{"haar"},
#'   \code{"d4"}, \code{"d6"}, \code{"d8"}, \code{"la16"}.
#' @param n_levels Integer scalar or NULL. MODWT levels. If NULL, chosen from CWT significance
#'   (when available) and constrained by stability and max-period caps.
#' @param boundary Character. Boundary handling method. Only \code{"periodic"} is supported.
#' @param include_smooth Logical. If TRUE, includes the smooth component SJ.
#' @param max_period_frac Numeric scalar in (0, 1]. No structure beyond this fraction of n.
#'   Default is \code{1/3}.
#' @param cwt_mode Character. CWT routine mode. One of \code{"fast"} or \code{"complete"}.
#' @param diagnostics Logical. If TRUE, attaches additional diagnostics (caps and covariance).
#'
#' @return A list with class \code{"wavelet_additive"}. The object contains:
#' - \code{cwt}: output from \code{analyze_wavelet_spectrum()}
#' - \code{components}: additive MODWT MRA components (matrix)
#' - \code{component_names}: component column names
#' - \code{periods}: representative MODWT periods per component
#' - \code{variance}, \code{variance_fraction}: component variance and fraction of total
#' - \code{n_levels}, \code{filter}, \code{filter_length}, \code{boundary}, \code{include_smooth},
#'   \code{max_period_frac}
#' - \code{has_significance}, \code{cwt_signif_periods}: CWT significance flags and periods
#' - \code{cwt_to_modwt_map}: MODWT levels implied by significant CWT periods (under caps)
#' - variance-accounting fields copied from \code{modwt_mra()}
#' - optional \code{diagnostics} list when requested
#'
#' @export
analyze_wavelet_additive <- function(
    series,
    signif = 0.90,
    noise = "red",
    min_period = 2,
    detrend = FALSE,
    filter = c("la8", "haar", "d4", "d6", "d8", "la16"),
    n_levels = NULL,
    boundary = "periodic",
    include_smooth = TRUE,
    max_period_frac = 1/3,
    cwt_mode = c("fast", "complete"),
    diagnostics = FALSE
) {
  series <- .as_numeric_series(series, "series", min_n = 16L)
  boundary <- .check_boundary(boundary)
  filter <- match.arg(filter)
  cwt_mode <- match.arg(cwt_mode)

  n <- length(series)
  filter_length <- .get_filter_length(filter)

  # --- CWT (diagnostic only) ---
  cwt_result <- analyze_wavelet_spectrum(
    series = series,
    signif = signif,
    noise = noise,
    min_period = min_period,
    detrend = detrend,
    mode = cwt_mode,
    diagnostics = diagnostics
  )

  # --- Determine allowable MODWT levels ---
  max_levels <- .recommended_max_levels(
    n = n,
    filter_length = filter_length,
    max_period_frac = max_period_frac
  )

  # --- Choose n_levels (optionally guided by CWT significance) ---
  if (is.null(n_levels)) {
    n_levels <- max(1L, max_levels - 1L)

    if (isTRUE(cwt_result$has_significance)) {
      idx <- cwt_result$signif_periods
      if (is.logical(idx)) idx <- which(idx)

      sig_periods <- cwt_result$period[idx]
      max_period_allowed <- floor(n * max_period_frac)

      sig_periods <- sig_periods[
        is.finite(sig_periods) & sig_periods >= 0 & sig_periods <= max_period_allowed
      ]

      if (length(sig_periods) > 0) {
        max_sig_period <- max(sig_periods)
        j_needed <- max(1L, floor(log2(max_sig_period)))
        n_levels <- min(max_levels, max(2L, as.integer(j_needed)))
      }
    }
  }

  n_levels <- .select_n_levels(n_levels, max_levels, default_backoff = 0L)

  # --- MODWT MRA (additive components) ---
  mra_result <- modwt_mra(
    x = series,
    filter = filter,
    n_levels = n_levels,
    boundary = boundary,
    include_smooth = include_smooth,
    max_period_frac = max_period_frac
  )

  # --- Map significant CWT periods to MODWT levels (under cap) ---
  cwt_to_modwt_map <- integer(0)
  cwt_signif_periods <- numeric(0)

  if (isTRUE(cwt_result$has_significance)) {
    idx <- cwt_result$signif_periods
    if (is.logical(idx)) idx <- which(idx)

    sig_periods <- cwt_result$period[idx]
    max_period_allowed <- floor(n * max_period_frac)

    sig_periods <- sig_periods[
      is.finite(sig_periods) & sig_periods >= 0 & sig_periods <= max_period_allowed
    ]

    cwt_signif_periods <- sig_periods

    if (length(sig_periods) > 0) {
      cwt_to_modwt_map <- vapply(sig_periods, function(p) {
        level <- floor(log2(p))
        level <- max(1L, min(level, mra_result$n_levels))
        as.integer(level)
      }, integer(1))

      cwt_to_modwt_map <- unique(cwt_to_modwt_map)
    }
  }

  comp_names <- colnames(mra_result$components)

  out <- list(
    cwt = cwt_result,
    components = mra_result$components,
    component_names = comp_names,
    periods = mra_result$periods,
    variance = mra_result$variance,
    variance_fraction = mra_result$variance_fraction,
    n_levels = mra_result$n_levels,
    filter = mra_result$filter,
    filter_length = mra_result$filter_length,
    boundary = mra_result$boundary,
    include_smooth = mra_result$include_smooth,
    max_period_frac = mra_result$max_period_frac,
    has_significance = isTRUE(cwt_result$has_significance),
    cwt_signif_periods = cwt_signif_periods,
    cwt_to_modwt_map = cwt_to_modwt_map,
    is_additive = TRUE,
    reconstruction_error = mra_result$reconstruction_error,
    cross_covariance_sum = mra_result$cross_covariance_sum,
    variance_identity_error = mra_result$variance_identity_error,
    total_variance = mra_result$total_variance,
    sum_component_variance = mra_result$sum_component_variance
  )

  if (isTRUE(diagnostics)) {
    out$diagnostics <- list(
      n_original = n,
      max_levels_final = max_levels,
      component_sd = sqrt(mra_result$variance),
      cov_matrix = stats::cov(mra_result$components)
    )
  }

  class(out) <- c("wavelet_additive", "list")
  out
}

#' Print method for additive wavelet analysis
#'
#' @description
#' Prints a compact summary of a \code{"wavelet_additive"} object.
#'
#' @param x Object of class \code{"wavelet_additive"}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.wavelet_additive <- function(x, ...) {
  cat("Additive Wavelet Analysis (CWT + MODWT MRA Hybrid)\n")
  cat(sprintf("  Series length:   %d\n", nrow(x$components)))
  cat(sprintf("  MODWT levels:    %d (%s filter, L=%d)\n", x$n_levels, x$filter, x$filter_length))
  cat(sprintf("  Components:      %d\n", ncol(x$components)))
  cat(sprintf("  max_period_frac: %.6f\n", x$max_period_frac))
  cat("\n")

  cat("Variance Decomposition:\n")

  period_labels <- sprintf("%.1f", x$periods)
  n_comp <- length(x$component_names)

  if (isTRUE(x$include_smooth) && x$component_names[n_comp] == paste0("S", x$n_levels)) {
    period_labels[n_comp] <- paste0(">", period_labels[n_comp])
  }

  var_df <- data.frame(
    Component = x$component_names,
    Period = period_labels,
    Variance = sprintf("%.6f", x$variance),
    Fraction = sprintf("%.1f%%", x$variance_fraction * 100)
  )
  print(var_df, row.names = FALSE)

  cat(sprintf("  Fractions sum to %.1f%% (deviation reflects cross-covariances)\n",
              sum(x$variance_fraction) * 100))
  cat("\n")

  cat("Variance Accounting:\n")
  cat(sprintf("  Var(x):                  %.8f\n", x$total_variance))
  cat(sprintf("  Sum var(components):      %.8f\n", x$sum_component_variance))
  cat(sprintf("  2 * sum cross-cov:        %.8f\n", 2 * x$cross_covariance_sum))
  cat(sprintf("  Identity error:           %.2e (target is near 0)\n", x$variance_identity_error))
  cat(sprintf("  Reconstruction error:     %.2e\n", x$reconstruction_error))
  cat("\n")

  if (isTRUE(x$has_significance)) {
    cat(sprintf("CWT Significance: %d significant period(s)\n", length(x$cwt_signif_periods)))
    if (length(x$cwt_signif_periods) > 0) {
      cat(sprintf("  Periods: %s\n", paste(round(x$cwt_signif_periods, 1), collapse = ", ")))
    }
    if (length(x$cwt_to_modwt_map) > 0) {
      cat(sprintf("  Maps to MODWT levels: %s\n", paste(x$cwt_to_modwt_map, collapse = ", ")))
    } else {
      cat("  Maps to MODWT levels: none under max_period_frac cap\n")
    }
  } else {
    cat("CWT Significance: none detected\n")
  }

  invisible(x)
}
