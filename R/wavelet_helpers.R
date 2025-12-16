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

  # Input validation
  if (!is.numeric(k) || !is.numeric(s) || !is.numeric(k0)) {
    stop("All inputs must be numeric")
  }

  if (length(s) != 1 || s <= 0) {
    stop("'s' must be a positive scalar")
  }

  if (length(k0) != 1 || k0 <= 0) {
    stop("'k0' must be a positive scalar")
  }

  nn <- length(k)

  # Only positive frequencies (one-sided spectrum)
  z <- as.numeric(k > 0)

  # Exponential term: Gaussian envelope in frequency
  expnt <- -((s * k - k0)^2 / 2) * z

  # Normalization factor (Torrence & Compo 1998, Eq. 4)
  norm <- sqrt(s * k[2]) * (pi^(-0.25)) * sqrt(nn)

  # Daughter wavelet in Fourier space
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
#' user-specified significant periods. Uses the inverse continuous wavelet
#' transform formula to convert wavelet coefficients back to time domain.
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
#' @param Cdelta Numeric scalar. Reconstruction constant for Morlet wavelet
#'   (default = 0.776). Depends on wavelet type; see Torrence & Compo (1998) Table 2.
#' @param w0_0 Numeric scalar. Normalization constant \eqn{\psi_0(0) = \pi^{-1/4}}
#'   for Morlet wavelet (default = \code{pi^(-1/4)}).
#'
#' @return A tibble with one column per significant period component.
#'   Column names are \code{Component_1}, \code{Component_2}, etc.
#'   Number of rows equals the length of the original time series.
#'
#' @details
#' The reconstruction formula (Torrence & Compo 1998, Eq. 11) is:
#' \deqn{x_n = \frac{\delta j \sqrt{\delta t}}{C_\delta \psi_0(0)}
#'       \sum_{j=0}^{J} \frac{\text{Re}(W_n(s_j))}{\sqrt{s_j}}}
#'
#' Each component represents the contribution of a specific frequency band
#' to the original time series. Components can be summed to approximate
#' the original signal (within the limits of the selected periods).
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
#'   Cdelta = 0.776
#' )
#'
#' # Components sum to approximate original signal
#' reconstructed <- rowSums(components)
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
                                       Cdelta = 0.776,
                                       w0_0 = pi^(-1/4)) {

  # Input validation
  if (!is.matrix(wave) && !is.array(wave)) {
    stop("'wave' must be a matrix or array")
  }

  if (!is.numeric(signif_periods) || length(signif_periods) == 0) {
    stop("'signif_periods' must be a non-empty numeric vector")
  }

  if (any(signif_periods < 1) || any(signif_periods > nrow(wave))) {
    stop("'signif_periods' indices must be between 1 and nrow(wave)")
  }

  if (!is.numeric(scale) || length(scale) != nrow(wave)) {
    stop("'scale' must be numeric with length equal to nrow(wave)")
  }

  if (dj <= 0 || dt <= 0 || variable_sd <= 0 || Cdelta <= 0 || w0_0 <= 0) {
    stop("Parameters dj, dt, variable_sd, Cdelta, and w0_0 must be positive")
  }

  num_periods <- length(signif_periods)
  n <- ncol(wave)

  # Preallocate component matrix
  COMPS <- matrix(0, nrow = n, ncol = num_periods)

  # Reconstruction factor (constant across all components)
  # Torrence & Compo (1998), Eq. 11
  recon_fac <- variable_sd * (dj * sqrt(dt) / (Cdelta * w0_0))

  # Reconstruct each component
  for (i in seq_len(num_periods)) {

    # Get scale index for current significant period
    cur_period_idx <- signif_periods[i]
    sj <- scale[cur_period_idx]

    # Extract real part of wavelet coefficients at this scale
    W <- Re(wave[cur_period_idx, , drop = FALSE])

    # Apply reconstruction formula: normalize by scale and multiply by factor
    W_scaled <- W / sqrt(sj)
    component <- as.numeric(recon_fac * W_scaled)

    COMPS[, i] <- component
  }

  # Convert to tibble with informative column names
  comps_tb <- tibble::as_tibble(
    COMPS,
    .name_repair = ~ paste0("Component_", seq_len(num_periods))
  )

  return(comps_tb)
}
