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


#' Get Wavelet Filter Length
#'
#' @description
#' Returns the filter length for a named wavelet filter.
#'
#' @param filter Character. Filter name.
#'
#' @return Integer. Filter length.
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

  if (!filter %in% names(lengths)) {
    stop("Unknown filter: ", filter, call. = FALSE)
  }

  lengths[[filter]]
}


#' Compute max allowed MODWT levels given a max-period fraction
#'
#' @description
#' Ensures the smooth threshold period 2^(J+1) does not exceed floor(n * max_period_frac).
#' This prevents representing variability at periods longer than a chosen fraction
#' of record length.
#'
#' @param n Integer. Series length.
#' @param max_period_frac Numeric in (0, 1]. Maximum allowed period as a fraction of n.
#'
#' @return Integer. Maximum allowed J (n_levels), at least 1.
#'
#' @keywords internal
.max_levels_by_period_frac <- function(n, max_period_frac) {
  if (!is.numeric(max_period_frac) || length(max_period_frac) != 1L ||
      !is.finite(max_period_frac) || max_period_frac <= 0 || max_period_frac > 1) {
    stop("'max_period_frac' must be in (0, 1].", call. = FALSE)
  }

  max_period_allowed <- floor(as.numeric(n) * max_period_frac)

  # Need at least 4 to allow J>=1 because 2^(J+1) with J=1 => 4
  if (!is.finite(max_period_allowed) || max_period_allowed < 4L) return(1L)

  Jmax <- floor(log2(max_period_allowed)) - 1L
  max(1L, as.integer(Jmax))
}


#' Maximal Overlap Discrete Wavelet Transform (MODWT)
#'
#' @description
#' Performs the Maximal Overlap Discrete Wavelet Transform using a specified
#' wavelet filter. MODWT is shift-invariant and provides additive multi-resolution
#' decomposition at dyadic scales (powers of 2).
#'
#' This function wraps \code{waveslim::modwt()} with additional input validation
#' and a standardized return structure.
#'
#' @param x Numeric vector. Input time series.
#' @param filter Character. Wavelet filter name. One of:
#'   \code{"haar"}, \code{"d4"}, \code{"d6"}, \code{"d8"},
#'   \code{"la8"}, \code{"la16"} (default: \code{"la8"}).
#' @param n_levels Integer. Number of decomposition levels. If \code{NULL},
#'   uses a conservative default based on series length and filter width.
#'   This is a stability heuristic, not a mathematical requirement of MODWT.
#' @param boundary Character. Boundary handling method. Currently only
#'   \code{"periodic"} is supported (circular convolution).
#'
#' @return A list of class \code{"modwt_result"} containing:
#' \describe{
#'   \item{coefficients}{List from waveslim::modwt with W1...WJ and VJ (scaling).}
#'   \item{filter}{Character. Filter name used.}
#'   \item{n_levels}{Integer. Number of decomposition levels (J).}
#'   \item{n}{Integer. Original series length.}
#'   \item{boundary}{Character. Boundary method used.}
#'   \item{filter_length}{Integer. Length of the wavelet filter.}
#' }
#'
#' @details
#' The number of decomposition levels is capped using a stability heuristic:
#' \code{min(floor(log2(n)), floor(log2((n-1)/(L-1)+1)))} where L is the filter
#' length. This ensures the effective filter support does not dominate the
#' series length, reducing boundary artifacts.
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

  if (!is.numeric(x)) stop("'x' must be a numeric vector", call. = FALSE)
  if (anyNA(x)) stop("'x' contains NA values", call. = FALSE)

  x <- as.numeric(x)
  n <- length(x)

  if (n < 8) stop("'x' must have at least 8 observations", call. = FALSE)

  filter <- match.arg(filter)

  if (!identical(boundary, "periodic")) {
    stop("Only 'periodic' boundary is supported.", call. = FALSE)
  }

  filter_length <- .get_filter_length(filter)

  max_levels_length <- floor(log2(n))
  max_levels_filter <- floor(log2((n - 1) / (filter_length - 1) + 1))
  max_levels <- max(1L, min(max_levels_length, max_levels_filter))

  if (is.null(n_levels)) {
    n_levels <- max(1L, max_levels - 1L)
  }

  if (!is.numeric(n_levels) || length(n_levels) != 1L || !is.finite(n_levels) || n_levels < 1) {
    stop("'n_levels' must be a positive integer", call. = FALSE)
  }

  n_levels <- as.integer(n_levels)

  if (n_levels > max_levels) {
    warning(
      sprintf(
        "Requested n_levels=%d exceeds recommended maximum=%d for n=%d with %s filter (L=%d). Using %d levels.",
        n_levels, max_levels, n, filter, filter_length, max_levels
      ),
      call. = FALSE
    )
    n_levels <- max_levels
  }

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
#' Reconstructs a time series from its MODWT coefficients.
#'
#' @param wt A list of class \code{"modwt_result"} from \code{\link{modwt_decompose}},
#'   or a raw waveslim modwt object.
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
    has_W <- any(grepl("^W[0-9]+$", names(wt)))
    has_V <- any(grepl("^V[0-9]+$", names(wt))) || ("VJ" %in% names(wt))

    if (!has_W || !has_V) {
      stop(
        "'wt' must be a 'modwt_result' object or valid waveslim modwt output ",
        "(requires W* and VJ/V<level> coefficients).",
        call. = FALSE
      )
    }
    coeffs <- wt

  } else {
    stop("'wt' must be a 'modwt_result' object or waveslim modwt output", call. = FALSE)
  }

  as.numeric(waveslim::imodwt(coeffs))
}


#' MODWT Multiresolution Analysis (MRA)
#'
#' @description
#' Extracts additive time-domain components from MODWT decomposition using
#' multiresolution analysis. Each component captures variance at a specific
#' dyadic scale band. The components sum exactly to the original series.
#'
#' @param x Numeric vector. Input time series.
#' @param filter Character. Wavelet filter name.
#' @param n_levels Integer. Number of decomposition levels (J). If \code{NULL},
#'   automatically determined based on series length, filter length, and max-period cap.
#' @param boundary Character. Boundary handling. Only \code{"periodic"} is supported.
#' @param include_smooth Logical. If \code{TRUE}, includes the smooth (approximation)
#'   component at the coarsest level (SJ).
#' @param max_period_frac Numeric in (0, 1]. Maximum allowed represented period as a
#'   fraction of record length. Default is 1/3 (no structure beyond n/3).
#'
#' @return A list with additive components, periods, variance accounting, and diagnostics.
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

  if (!is.numeric(x)) stop("'x' must be a numeric vector", call. = FALSE)
  if (anyNA(x)) stop("'x' contains NA values", call. = FALSE)

  x <- as.numeric(x)
  n <- length(x)

  if (n < 8) stop("'x' must have at least 8 observations", call. = FALSE)

  filter <- match.arg(filter)

  if (!identical(boundary, "periodic")) {
    stop("Only 'periodic' boundary is supported.", call. = FALSE)
  }

  filter_length <- .get_filter_length(filter)

  # Base stability heuristic
  max_levels_length <- floor(log2(n))
  max_levels_filter <- floor(log2((n - 1) / (filter_length - 1) + 1))
  max_levels <- max(1L, min(max_levels_length, max_levels_filter))

  # Additional cap: no structure beyond a fraction of record length
  max_levels_period <- .max_levels_by_period_frac(n, max_period_frac)
  max_levels <- min(max_levels, max_levels_period)

  if (is.null(n_levels)) {
    n_levels <- max(1L, max_levels - 1L)
  }

  if (!is.numeric(n_levels) || length(n_levels) != 1L || !is.finite(n_levels) || n_levels < 1) {
    stop("'n_levels' must be a positive integer", call. = FALSE)
  }

  n_levels <- as.integer(min(n_levels, max_levels))
  n_levels <- max(1L, n_levels)

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

  # Approximate Fourier periods:
  # Detail j ~ band [2^j, 2^(j+1)] -> report geometric mean 2^j * sqrt(2)
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


#' Additive Wavelet Analysis with CWT Diagnostics
#'
#' @description
#' Combines CWT analysis (visualization + significance) with MODWT MRA (additive components).
#' CWT is used only for diagnostics; MODWT MRA provides additive component extraction.
#'
#' @param series Numeric vector. Input series (regularly spaced, no missing values).
#' @param signif Numeric scalar in (0, 1). Significance level for CWT.
#' @param noise Character. Background noise model for CWT significance testing ("white" or "red").
#' @param min_period Numeric. Minimum Fourier period for CWT.
#' @param detrend Logical. If TRUE, removes linear trend before CWT.
#' @param filter Character. MODWT filter.
#' @param n_levels Integer or NULL. MODWT levels. If NULL, selected from CWT + caps.
#' @param boundary Character. Only "periodic" supported.
#' @param include_smooth Logical. If TRUE, include smooth component SJ.
#' @param max_period_frac Numeric in (0, 1]. No structure beyond this fraction of n (default 1/3).
#' @param cwt_mode Character. "fast" or "complete" for CWT routine.
#' @param diagnostics Logical. If TRUE, attach extra diagnostics (cov matrix, caps, etc.).
#'
#' @return A list of class "wavelet_additive".
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

  if (!is.numeric(series)) stop("'series' must be numeric", call. = FALSE)
  if (anyNA(series)) stop("'series' contains missing values", call. = FALSE)
  if (length(series) < 16) stop("'series' must have at least 16 observations", call. = FALSE)

  filter <- match.arg(filter)
  cwt_mode <- match.arg(cwt_mode)

  if (!identical(boundary, "periodic")) {
    stop("Only 'periodic' boundary is supported.", call. = FALSE)
  }

  n <- length(series)
  series <- as.numeric(series)

  # --- CWT Analysis (diagnostic only) ---
  # NOTE: analyze_wavelet_spectrum() must exist elsewhere in your package.
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
  filter_length <- .get_filter_length(filter)

  max_levels_length <- floor(log2(n))
  max_levels_filter <- floor(log2((n - 1) / (filter_length - 1) + 1))
  max_levels <- max(1L, min(max_levels_length, max_levels_filter))

  # Additional cap: no structure beyond a fraction of record length
  max_levels_period <- .max_levels_by_period_frac(n, max_period_frac)
  max_levels <- min(max_levels, max_levels_period)

  if (is.null(n_levels)) {
    if (isTRUE(cwt_result$has_significance)) {
      idx <- cwt_result$signif_periods
      if (is.logical(idx)) idx <- which(idx)
      sig_periods <- cwt_result$period[idx]

      # Optional: ignore CWT periods beyond the allowed cap
      max_period_allowed <- floor(n * max_period_frac)
      sig_periods <- sig_periods[is.finite(sig_periods) & sig_periods <= max_period_allowed]

      if (length(sig_periods) > 0) {
        max_sig_period <- max(sig_periods)

        # Need level j where 2^j <= p < 2^(j+1), so j = floor(log2(p))
        j_needed <- max(1L, floor(log2(max_sig_period)))
        n_levels <- min(max_levels, max(2L, j_needed))
      } else {
        # No usable significant periods under cap -> default conservative
        n_levels <- max(1L, max_levels - 1L)
      }

    } else {
      n_levels <- max(1L, max_levels - 1L)
    }
  }

  n_levels <- as.integer(min(n_levels, max_levels))
  n_levels <- max(1L, n_levels)

  # --- MODWT MRA (additive components) ---
  mra_result <- modwt_mra(
    x = series,
    filter = filter,
    n_levels = n_levels,
    boundary = boundary,
    include_smooth = include_smooth,
    max_period_frac = max_period_frac
  )

  # --- Map CWT significant periods to MODWT levels (under cap) ---
  cwt_to_modwt_map <- integer(0)

  if (isTRUE(cwt_result$has_significance)) {
    idx <- cwt_result$signif_periods
    if (is.logical(idx)) idx <- which(idx)
    sig_periods <- cwt_result$period[idx]

    max_period_allowed <- floor(n * max_period_frac)
    sig_periods <- sig_periods[is.finite(sig_periods) & sig_periods <= max_period_allowed]

    if (length(sig_periods) > 0) {
      cwt_to_modwt_map <- vapply(sig_periods, function(p) {
        level <- floor(log2(p))
        level <- max(1L, min(level, n_levels))
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
    cwt_signif_periods = if (isTRUE(cwt_result$has_significance)) {
      idx <- cwt_result$signif_periods
      if (is.logical(idx)) idx <- which(idx)
      cwt_result$period[idx]
    } else {
      numeric(0)
    },
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
      max_levels_length = max_levels_length,
      max_levels_filter = max_levels_filter,
      max_levels_period = max_levels_period,
      max_levels_final = max_levels,
      component_sd = sqrt(mra_result$variance),
      cov_matrix = stats::cov(mra_result$components)
    )
  }

  class(out) <- c("wavelet_additive", "list")
  out
}


#' Print Method for Additive Wavelet Analysis
#'
#' @param x Object of class \code{"wavelet_additive"}.
#' @param ... Ignored.
#' @return Invisibly returns x.
#'
#' @export
print.wavelet_additive <- function(x, ...) {

  cat("Additive Wavelet Analysis (CWT + MODWT MRA Hybrid)\n")
  cat(sprintf("  Series length: %d\n", nrow(x$components)))
  cat(sprintf("  MODWT levels:  %d (%s filter, L=%d)\n",
              x$n_levels, x$filter, x$filter_length))
  cat(sprintf("  Components:    %d\n", ncol(x$components)))
  cat(sprintf("  max_period_frac: %.4f\n", x$max_period_frac))
  cat("\n")

  cat("Variance Decomposition:\n")

  period_labels <- sprintf("%.1f", x$periods)
  n_comp <- length(x$component_names)
  if (x$include_smooth && x$component_names[n_comp] == paste0("S", x$n_levels)) {
    period_labels[n_comp] <- paste0(">", period_labels[n_comp])
  }

  var_df <- data.frame(
    Component = x$component_names,
    Period = period_labels,
    Variance = sprintf("%.6f", x$variance),
    Fraction = sprintf("%.1f%%", x$variance_fraction * 100)
  )
  print(var_df, row.names = FALSE)
  cat(sprintf("  (Fractions sum to %.1f%% due to cross-covariances)\n",
              sum(x$variance_fraction) * 100))
  cat("\n")

  cat("Variance Accounting:\n")
  cat(sprintf("  Original Var(x):        %.8f\n", x$total_variance))
  cat(sprintf("  Sum of component var:   %.8f\n", x$sum_component_variance))
  cat(sprintf("  2 * cross-cov sum:      %.8f\n", 2 * x$cross_covariance_sum))
  cat(sprintf("  Identity error:         %.2e (should be ~0)\n", x$variance_identity_error))
  cat(sprintf("  Reconstruction error:   %.2e\n", x$reconstruction_error))
  cat("\n")

  if (isTRUE(x$has_significance)) {
    cat(sprintf("CWT Significance: %d significant period(s) detected\n",
                length(x$cwt_signif_periods)))
    cat(sprintf("  Periods: %s\n", paste(round(x$cwt_signif_periods, 1), collapse = ", ")))
    if (length(x$cwt_to_modwt_map) > 0) {
      cat(sprintf("  Maps to MODWT levels: %s\n", paste(x$cwt_to_modwt_map, collapse = ", ")))
    } else {
      cat("  Maps to MODWT levels: (none under max_period_frac cap)\n")
    }
  } else {
    cat("CWT Significance: No significant periods detected\n")
  }

  invisible(x)
}
