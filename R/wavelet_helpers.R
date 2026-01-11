# ==============================================================================
# Wavelet Utility Functions
# ==============================================================================
# This file contains utility functions for wavelet analysis and filtering.
# New additions for filter_warm_simulations support:
# - fill_nearest()
# - extract_signif_curve()
# - gws_regrid()
# ==============================================================================

#' Fill NA values with nearest non-NA value
#'
#' @description
#' Forward and backward fills NA values in a vector using the nearest
#' non-NA value. Useful for handling edge effects in interpolated data.
#'
#' @param v Numeric vector possibly containing NAs
#'
#' @return Numeric vector with NAs filled
#'
#' @keywords internal
#' @export
fill_nearest <- function(v) {
  if (all(is.na(v))) return(v)

  # Forward fill from first non-NA
  if (is.na(v[1])) {
    first <- which(!is.na(v))[1]
    if (!is.na(first)) v[1:(first - 1)] <- v[first]
  }

  # Forward fill
  for (i in 2:length(v)) {
    if (is.na(v[i])) v[i] <- v[i - 1]
  }

  # Backward fill
  for (i in (length(v) - 1):1) {
    if (is.na(v[i])) v[i] <- v[i + 1]
  }

  v
}

#' Extract significance curve from wavelet analysis result
#'
#' @description
#' Searches for GWS significance curve in wavelet analysis output.
#' Tries multiple possible field names.
#'
#' @param wv List output from wavelet_spectral_analysis()
#'
#' @return Numeric vector of significance values, or NULL if not found
#'
#' @keywords internal
#' @export
extract_signif_curve <- function(wv) {
  # Try multiple possible names
  cand <- c("gws_signif", "gws_sig", "signif_gws", "signif_gws_curve",
            "gws_significance", "signif")

  for (nm in cand) {
    if (!is.null(wv[[nm]]) && is.numeric(wv[[nm]])) {
      return(as.numeric(wv[[nm]]))
    }
  }

  NULL
}

#' Regrid GWS to target period vector
#'
#' @description
#' Interpolates GWS from wavelet analysis onto a target period grid.
#' Uses linear interpolation and fills edge NAs.
#'
#' @param wv List output from wavelet_spectral_analysis()
#' @param target_period Numeric vector of target periods
#' @param use_unmasked Logical. If TRUE and available, use gws_unmasked
#'
#' @return Numeric vector of GWS values on target_period grid
#'
#' @keywords internal
#' @export
gws_regrid <- function(wv, target_period, use_unmasked = FALSE) {
  # Get GWS
  g <- if (isTRUE(use_unmasked) && !is.null(wv$gws_unmasked)) {
    wv$gws_unmasked
  } else {
    wv$gws
  }

  # Get period
  p <- wv$gws_period

  # Validate
  if (is.null(g) || !is.numeric(g)) {
    stop("wavelet_spectral_analysis() missing numeric $gws.", call. = FALSE)
  }
  if (is.null(p) || !is.numeric(p)) {
    stop("wavelet_spectral_analysis() missing numeric $gws_period.", call. = FALSE)
  }

  g <- as.numeric(g)
  p <- as.numeric(p)

  # Remove non-finite values
  ok <- is.finite(p) & is.finite(g)
  p <- p[ok]
  g <- g[ok]

  # Handle edge case: too few points
  if (length(g) < 2L) {
    out <- rep(if (length(g) == 1L) g[1] else NA_real_, length(target_period))
    return(fill_nearest(out))
  }

  # Interpolate
  out <- stats::approx(x = p, y = g, xout = target_period, rule = 2)$y
  fill_nearest(out)
}

# ==============================================================================
# Original wavelet_helpers.R functions below
# ==============================================================================

#' Morlet Wavelet in Fourier Domain
#'
#' @description
#' Computes the Morlet wavelet function in Fourier space for a given scale.
#' This is the "daughter" wavelet used in continuous wavelet transform analysis.
#' The function applies normalization following Torrence & Compo (1998).
#'
#' @param k Numeric vector. Wave number vector (angular frequency).
#' @param s Numeric scalar. Scale parameter (inverse of frequency).
#' @param k0 Numeric scalar. Omega0 parameter for the Morlet wavelet (default = 6).
#'   This value balances time and frequency localization.
#'
#' @return Numeric vector (complex). The Morlet wavelet in Fourier space,
#'   with length equal to length(k).
#'
#' @details
#' The Morlet wavelet is defined in Fourier space as:
#' psi(s k) = pi^(-1/4) * sqrt(s k_2 n) * exp(-(s k - k_0)^2 / 2)
#' where the exponential term is applied only to positive frequencies.
#'
#' The normalization factor ensures energy conservation and proper inverse
#' transform reconstruction.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' Bulletin of the American Meteorological Society, 79(1), 61-78.
#'
#' @keywords internal
#' @export
morlet_wavelet <- function(k, s, k0 = 6) {

  # Basic validation
  if (!is.numeric(k) || !is.numeric(s) || !is.numeric(k0)) stop("All inputs must be numeric")
  if (length(s) != 1 || s <= 0) stop("'s' must be a positive scalar")
  if (length(k0) != 1 || k0 <= 0) stop("'k0' must be a positive scalar")
  if (length(k) < 2) stop("'k' must have length >= 2")

  if (abs(k[1]) > 1e-10) stop("'k' must start at zero frequency (standard FFT ordering)")
  if (k[2] <= 0) stop("'k[2]' must be positive (frequency spacing)")
  if (any(!is.finite(k))) stop("'k' contains non-finite values (NA, NaN, or Inf)")

  nn <- length(k)
  z <- as.numeric(k > 0)
  expnt <- -((s * k - k0)^2 / 2) * z

  # TC98 Eq. 4 normalization (TWO-SIDED / standard)
  norm <- sqrt(s * k[2]) * (pi^(-0.25)) * sqrt(nn)

  daughter <- norm * exp(expnt) * z
  daughter
}

#' Compute Morlet Wavelet Parameters
#'
#' @description
#' Calculates key parameters for the Morlet wavelet transform, including
#' the Fourier factor (scale-to-period conversion), cone of influence (COI),
#' and minimum degrees of freedom.
#'
#' @param k0 Numeric scalar. Omega0 parameter for the Morlet wavelet (default = 6).
#'   Standard choice is 6, which gives ~6 oscillations per wavelet.
#'
#' @return Named numeric vector with three elements:
#'   fourier_factor: Conversion factor from wavelet scale to Fourier period
#'   coi: Cone of influence e-folding time
#'   dofmin: Minimum degrees of freedom (2 for Morlet)
#'
#' @details
#' The Fourier factor allows conversion between wavelet scale s and
#' equivalent Fourier period T:
#' T = (4*pi / (k_0 + sqrt(2 + k_0^2))) * s
#'
#' The cone of influence (COI) defines the e-folding time for edge effects:
#' COI = fourier_factor / sqrt(2)
#'
#' Values outside the COI are subject to edge artifacts and should be
#' interpreted with caution.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' Bulletin of the American Meteorological Society, 79(1), 61-78.
#' See Table 1 and Section 3f.
#'
#' @examples
#' # Standard Morlet wavelet parameters
#' params <- morlet_parameters(k0 = 6)
#' print(params)
#' # fourier_factor: 1.03
#' # coi: 0.73
#' # dofmin: 2
#'
#' @keywords internal
#' @export
morlet_parameters <- function(k0 = 6) {

  # Input validation
  if (!is.numeric(k0) || length(k0) != 1 || k0 <= 0) {
    stop("'k0' must be a positive numeric scalar")
  }

  # Fourier factor: converts scale to period (TC98 Table 1)
  fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0^2))

  # Cone of influence e-folding time (TC98 Table 1)
  coi <- fourier_factor / sqrt(2)

  # Degrees of freedom for Morlet wavelet (TC98 Table 1)
  dofmin <- 2

  return(c(
    fourier_factor = fourier_factor,
    coi = coi,
    dofmin = dofmin
  ))
}

#' Extract and Reconstruct Wavelet Components from Significant Periods
#'
#' @description
#' Reconstructs time series components from the wavelet transform for
#' user-specified significant periods. Optionally includes a residual component
#' representing all non-significant scales, ensuring complete signal decomposition.
#'
#' @param wave Complex matrix. Wavelet transform coefficients with dimensions
#'   (scales x time). Output from continuous wavelet transform.
#' @param signif_periods Integer vector. Indices of significant period scales
#'   to reconstruct. Each index corresponds to a row in wave.
#' @param scale Numeric vector. Scale values corresponding to rows of wave.
#' @param dj Numeric scalar. Scale resolution parameter used in wavelet transform
#'   (default = 0.25). Smaller values give finer resolution.
#' @param dt Numeric scalar. Time step of the original series (default = 1).
#' @param variable_sd Numeric scalar. Standard deviation of the original time series
#'   before standardization. Used to restore original units.
#' @param variable_mean Numeric scalar. Mean of the original time series
#'   before standardization (default = 0). Added back to noise component to ensure
#'   complete reconstruction: rowSums(output) = original_signal.
#' @param Cdelta Numeric scalar. Reconstruction constant for Morlet wavelet
#'   (default = 0.776). Depends on wavelet type; see Torrence & Compo (1998) Table 2.
#' @param w0_0 Numeric scalar. Normalization constant psi_0(0) = pi^(-1/4)
#'   for Morlet wavelet (default = pi^(-1/4)).
#' @param include_residual Logical. If TRUE (default), includes a "Noise" component
#'   representing the sum of all non-significant scales, ensuring complete
#'   signal decomposition where rowSums(output) equals the reconstructed signal.
#'
#' @return A matrix with columns:
#'   Component_1, Component_2, ...: One column per significant period
#'   Noise: (if include_residual = TRUE) Sum of all non-significant scales
#'   Number of rows equals the length of the original time series.
#'   Property: rowSums(output) approximately equals the original signal.
#'
#' @details
#' The reconstruction formula (Torrence & Compo 1998, Eq. 11) is:
#'
#' x_n = (delta_j * sqrt(delta_t)) / (C_delta * psi_0(0)) *
#'       sum over j of Re(W_n(s_j)) / sqrt(s_j)
#'
#' This function performs a complete wavelet decomposition by:
#' 1. Extracting each significant period as a separate component
#' 2. Summing all remaining (non-significant) scales into a "Noise" component
#' 3. Ensuring: original_signal = sum(significant_components) + Noise
#'
#' This differs from partial reconstruction, which would only sum significant
#' scales and lose information from non-significant scales.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' Bulletin of the American Meteorological Society, 79(1), 61-78.
#' See Section 3i and Table 2.
#'
#' @examples
#' \dontrun{
#' # After performing wavelet analysis
#' components <- extract_wavelet_components(
#'   wave = wave_transform,
#'   signif_periods = c(5, 10, 20),  # indices of significant scales
#'   scale = scale_vector,
#'   dj = 0.25,
#'   dt = 1,
#'   variable_sd = sd(original_data),
#'   variable_mean = mean(original_data),
#'   Cdelta = 0.776,
#'   include_residual = TRUE
#' )
#'
#' # Verify complete reconstruction
#' reconstructed <- rowSums(components)
#' plot(original_data)
#' lines(reconstructed, col = "red")
#'
#' # Components include:
#' # - Component_1: First significant period
#' # - Component_2: Second significant period
#' # - Component_3: Third significant period
#' # - Noise: All non-significant scales combined + mean
#' }
#'
#' @keywords internal
#' @export
extract_wavelet_components <- function(wave,
                                       signif_periods,
                                       scale,
                                       dj = 0.25,
                                       dt = 1,
                                       variable_sd,
                                       variable_mean = 0,
                                       Cdelta = 0.776,
                                       w0_0 = pi^(-1/4),
                                       include_residual = TRUE) {

  if (!is.matrix(wave) && !is.array(wave)) stop("'wave' must be a matrix or array")
  if (!is.numeric(signif_periods) || length(signif_periods) == 0) stop("'signif_periods' must be a non-empty numeric vector")
  if (!all(signif_periods == as.integer(signif_periods))) stop("'signif_periods' must contain integer indices")
  if (any(signif_periods < 1) || any(signif_periods > nrow(wave))) stop("'signif_periods' indices must be between 1 and nrow(wave)")
  if (anyDuplicated(signif_periods)) {
    warning("'signif_periods' contains duplicate indices; duplicates will be removed")
    signif_periods <- unique(signif_periods)
  }

  if (!is.numeric(scale) || length(scale) != nrow(wave)) stop("'scale' must be numeric with length equal to nrow(wave)")
  if (!is.numeric(dj) || length(dj) != 1 || dj <= 0) stop("'dj' must be a positive scalar")
  if (!is.numeric(dt) || length(dt) != 1 || dt <= 0) stop("'dt' must be a positive scalar")
  if (!is.numeric(variable_sd) || length(variable_sd) != 1 || variable_sd <= 0) stop("'variable_sd' must be a positive scalar")
  if (!is.numeric(variable_mean) || length(variable_mean) != 1) stop("'variable_mean' must be a numeric scalar")
  if (!is.numeric(Cdelta) || length(Cdelta) != 1 || Cdelta <= 0) stop("'Cdelta' must be a positive scalar")
  if (!is.numeric(w0_0) || length(w0_0) != 1 || w0_0 <= 0) stop("'w0_0' must be a positive scalar")
  if (!is.logical(include_residual) || length(include_residual) != 1) stop("'include_residual' must be a single logical value")

  num_periods <- length(signif_periods)
  n_time <- ncol(wave)
  n_scales <- nrow(wave)

  n_components <- num_periods + (if (include_residual) 1L else 0L)
  COMPS <- matrix(0, nrow = n_time, ncol = n_components)

  # TC98 Eq. 11 reconstruction factor (TWO-SIDED / standard)
  recon_fac <- variable_sd * (dj * sqrt(dt) / (Cdelta * w0_0))

  inv_sqrt_scale <- 1 / sqrt(scale)
  Ww_all <- Re(wave) * inv_sqrt_scale

  for (i in seq_len(num_periods)) {
    j <- signif_periods[i]
    COMPS[, i] <- recon_fac * Ww_all[j, ]
  }

  if (include_residual) {
    nonsig_idx <- setdiff(seq_len(n_scales), signif_periods)

    if (length(nonsig_idx) > 0) {
      COMPS[, num_periods + 1L] <- recon_fac * colSums(Ww_all[nonsig_idx, , drop = FALSE])
    } else {
      COMPS[, num_periods + 1L] <- 0
      warning("All scales are marked as significant; residual component will be zero.")
    }

    COMPS[, num_periods + 1L] <- COMPS[, num_periods + 1L] + variable_mean
  }

  colnames(COMPS) <- c(
    paste0("Component_", seq_len(num_periods)),
    if (include_residual) "Noise" else NULL
  )

  attr(COMPS, "signif_periods") <- signif_periods
  attr(COMPS, "n_scales") <- n_scales
  attr(COMPS, "n_significant") <- num_periods
  attr(COMPS, "reconstruction_complete") <- include_residual

  COMPS
}
