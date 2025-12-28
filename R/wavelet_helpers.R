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
#'   with length equal to \code{length(k)}.
#'
#' @details
#' The Morlet wavelet is defined in Fourier space as:
#' \deqn{\psi(s k) = \pi^{-1/4} \sqrt{s k_2 n} \exp(-\frac{(s k - k_0)^2}{2})}
#' where the exponential term is applied only to positive frequencies.
#'
#' The normalization factor ensures energy conservation and proper inverse
#' transform reconstruction.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
#'
#' @keywords internal
#' @export
morlet_wavelet <- function(k, s, k0 = 6) {

  # Basic validation
  if (!is.numeric(k) || !is.numeric(s) || !is.numeric(k0)) {
    stop("All inputs must be numeric")
  }

  if (length(s) != 1 || s <= 0) {
    stop("'s' must be a positive scalar")
  }

  if (length(k0) != 1 || k0 <= 0) {
    stop("'k0' must be a positive scalar")
  }

  if (length(k) < 2) {
    stop("'k' must have length >= 2")
  }

  # Simple critical checks for FFT frequency vector
  if (abs(k[1]) > 1e-10) {
    stop("'k' must start at zero frequency (standard FFT ordering)")
  }

  if (k[2] <= 0) {
    stop("'k[2]' must be positive (frequency spacing)")
  }

  # Check for non-finite values
  if (any(!is.finite(k))) {
    stop("'k' contains non-finite values (NA, NaN, or Inf)")
  }

  # Morlet wavelet computation
  nn <- length(k)
  z <- as.numeric(k > 0)
  expnt <- -((s * k - k0)^2 / 2) * z

  # Normalization following TC98 Eq. 4
  # Factor breakdown:
  #   pi^(-1/4): standard Morlet normalization
  #   sqrt(2 * s * dk): scale and frequency spacing
  #   sqrt(n): FFT normalization
  # NOTE: The sqrt(2) factor is critical for energy conservation
  norm <- sqrt(2 * s * k[2]) * (pi^(-0.25)) * sqrt(nn)

  daughter <- norm * exp(expnt) * z

  return(daughter)
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
#'   \describe{
#'     \item{fourier_factor}{Conversion factor from wavelet scale to Fourier period}
#'     \item{coi}{Cone of influence e-folding time}
#'     \item{dofmin}{Minimum degrees of freedom (2 for Morlet)}
#'   }
#'
#' @details
#' The Fourier factor allows conversion between wavelet scale \eqn{s} and
#' equivalent Fourier period \eqn{T}:
#' \deqn{T = \frac{4\pi}{k_0 + \sqrt{2 + k_0^2}} \cdot s}
#'
#' The cone of influence (COI) defines the e-folding time for edge effects:
#' \deqn{COI = \frac{\text{fourier\_factor}}{\sqrt{2}}}
#'
#' Values outside the COI are subject to edge artifacts and should be
#' interpreted with caution.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
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
#'   to reconstruct. Each index corresponds to a row in \code{wave}.
#' @param scale Numeric vector. Scale values corresponding to rows of \code{wave}.
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
#'   for Morlet wavelet (default = \code{pi^(-1/4)}).
#' @param include_residual Logical. If TRUE (default), includes a "Noise" component
#'   representing the sum of all non-significant scales, ensuring complete
#'   signal decomposition where rowSums(output) equals the reconstructed signal.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{Component_1, Component_2, ...}: One column per significant period
#'     \item \code{Noise}: (if include_residual = TRUE) Sum of all non-significant scales
#'   }
#'   Number of rows equals the length of the original time series.
#'   Property: \code{rowSums(output)} approximately equals the original signal.
#'
#' @details
#' The reconstruction formula (Torrence & Compo 1998, Eq. 11) is:
#'
#' x_n = (delta_j * sqrt(delta_t)) / (C_delta * psi_0(0)) *
#'       sum over j of Re(W_n(s_j)) / sqrt(s_j)
#'
#' This function performs a complete wavelet decomposition by:
#' \enumerate{
#'   \item Extracting each significant period as a separate component
#'   \item Summing all remaining (non-significant) scales into a "Noise" component
#'   \item Ensuring: original_signal = sum(significant_components) + Noise
#' }
#'
#' This differs from partial reconstruction, which would only sum significant
#' scales and lose information from non-significant scales.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
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
#' @importFrom tibble as_tibble
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

  # --- Validation ---
  if (!is.matrix(wave) && !is.array(wave)) stop("'wave' must be a matrix or array")

  if (!is.numeric(signif_periods) || length(signif_periods) == 0) {
    stop("'signif_periods' must be a non-empty numeric vector")
  }

  if (!all(signif_periods == as.integer(signif_periods))) {
    stop("'signif_periods' must contain integer indices")
  }

  if (any(signif_periods < 1) || any(signif_periods > nrow(wave))) {
    stop("'signif_periods' indices must be between 1 and nrow(wave)")
  }

  if (anyDuplicated(signif_periods)) {
    warning("'signif_periods' contains duplicate indices; duplicates will be removed")
    signif_periods <- unique(signif_periods)
  }

  if (!is.numeric(scale) || length(scale) != nrow(wave)) {
    stop("'scale' must be numeric with length equal to nrow(wave)")
  }

  if (!is.numeric(dj) || length(dj) != 1 || dj <= 0) stop("'dj' must be a positive scalar")
  if (!is.numeric(dt) || length(dt) != 1 || dt <= 0) stop("'dt' must be a positive scalar")

  if (!is.numeric(variable_sd) || length(variable_sd) != 1 || variable_sd <= 0) {
    stop("'variable_sd' must be a positive scalar")
  }

  if (!is.numeric(variable_mean) || length(variable_mean) != 1) {
    stop("'variable_mean' must be a numeric scalar")
  }

  if (!is.numeric(Cdelta) || length(Cdelta) != 1 || Cdelta <= 0) stop("'Cdelta' must be a positive scalar")
  if (!is.numeric(w0_0) || length(w0_0) != 1 || w0_0 <= 0) stop("'w0_0' must be a positive scalar")

  if (!is.logical(include_residual) || length(include_residual) != 1) {
    stop("'include_residual' must be a single logical value")
  }

  # --- Setup ---
  num_periods <- length(signif_periods)
  n_time <- ncol(wave)
  n_scales <- nrow(wave)

  n_components <- num_periods + (if (include_residual) 1L else 0L)
  COMPS <- matrix(0, nrow = n_time, ncol = n_components)

  # --- Reconstruction (TC98 Eq. 11) ---
  # One-sided convention compensation consistent with your morlet_wavelet() sqrt(2)
  Cdelta_effective <- Cdelta * sqrt(2)
  recon_fac <- variable_sd * (dj * sqrt(dt) / (Cdelta_effective * w0_0))

  inv_sqrt_scale <- 1 / sqrt(scale)
  Ww_all <- Re(wave) * inv_sqrt_scale  # scales x time

  # Significant components (single scale indices)
  for (i in seq_len(num_periods)) {
    j <- signif_periods[i]
    COMPS[, i] <- recon_fac * Ww_all[j, ]
  }

  # Optional residual (all non-significant scales)
  if (include_residual) {
    nonsig_idx <- setdiff(seq_len(n_scales), signif_periods)

    if (length(nonsig_idx) > 0) {
      COMPS[, num_periods + 1L] <- recon_fac * colSums(Ww_all[nonsig_idx, , drop = FALSE])
    } else {
      COMPS[, num_periods + 1L] <- 0
      warning("All scales are marked as significant; residual component will be zero.")
    }

    # Put mean into residual only if you want full-signal reconstruction from helper alone
    COMPS[, num_periods + 1L] <- COMPS[, num_periods + 1L] + variable_mean
  }

  # Names
  colnames(COMPS) <- c(
    paste0("Component_", seq_len(num_periods)),
    if (include_residual) "Noise" else NULL
  )

  # Metadata
  attr(COMPS, "signif_periods") <- signif_periods
  attr(COMPS, "n_scales") <- n_scales
  attr(COMPS, "n_significant") <- num_periods
  attr(COMPS, "reconstruction_complete") <- include_residual

  COMPS
}

