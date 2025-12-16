#' Wavelet Spectral Analysis with Significance Testing
#'
#' Performs continuous wavelet transform analysis using the Morlet mother wavelet
#' on a time series. Computes the global wavelet spectrum with significance levels
#' and optionally reconstructs signal components from significant frequency bands.
#'
#' The function tests wavelet power against a background noise spectrum (white or red)
#' and identifies periods with statistically significant power. Component extraction
#' reconstructs the time series contributions from significant frequency bands.
#'
#' @param variable Numeric vector representing the time series (e.g., annual precipitation).
#'   Must have at least 32 observations for reliable wavelet analysis.
#' @param signif.level Numeric between 0 and 1. Significance level for test (default = 0.90).
#' @param noise.type Character: "white" or "red" (default). Background noise model for
#'   significance testing. For "red", lag-1 autocorrelation is computed from the data.
#' @param period.lower.limit Numeric; minimum period (in time units) to consider
#'   significant (default = 2). Filters out very high-frequency noise.
#'
#' @return A list with:
#'   \item{GWS}{Numeric vector: Global Wavelet Spectrum (average power across time)}
#'   \item{GWS_signif}{Numeric vector: Significance threshold for GWS at each period}
#'   \item{GWS_period}{Numeric vector: Period (in time units) corresponding to each scale}
#'   \item{signif_periods}{Integer vector: Indices of periods with significant power, or NULL if none}
#'   \item{COMPS}{Data frame: Reconstructed components for each significant period band,
#'     plus residual NOISE component. NULL if no significant periods found.}
#'
#' @details
#' This implementation uses the Morlet mother wavelet with omega0 = 6. Significance
#' testing follows Torrence and Compo (1998). For red noise, the lag-1 autocorrelation
#' is estimated from the input time series. Cone of influence and edge effects are
#' accounted for during significance calculations.
#'
#' @references
#' Torrence C and Compo GP (1998): A Practical Guide to Wavelet Analysis.
#' Bull. Amer. Meteor. Soc., 79, 61-78.
#'
#' @examples
#' set.seed(123)
#' # Simulate AR(1) process with 8-year periodicity
#' years <- 1:128
#' signal <- 2 * sin(2 * pi * years / 8) + arima.sim(n = 128, model = list(ar = 0.7))
#' res <- wavelet_spectral_analysis(signal, noise.type = "red")
#'
#' # Display significant periods
#' if (!is.null(res$signif_periods)) {
#'   cat("Significant periods (years):",
#'       round(res$GWS_period[res$signif_periods], 2), "\n")
#' }
#'
#' @importFrom stats var sd cor fft simulate
#' @import tibble
#' @export
wavelet_spectral_analysis <- function(variable,
                            signif.level = 0.90,
                            noise.type = "red",
                            period.lower.limit = 2,
                            detrend = FALSE) {

  # --- Input Validation ---

  # Check variable is numeric and sufficient length
  if (!is.numeric(variable)) {
    stop("Variable must be numeric.")
  }

  if (anyNA(variable)) {
    stop("Variable contains missing values. Please remove or interpolate NAs.")
  }

  # Check noise type
  if (!(noise.type %in% c("white", "red"))) {
    stop("noise.type must be 'white' or 'red'.")
  }

  # Check significance level
  if (!is.numeric(signif.level) || signif.level <= 0 || signif.level >= 1) {
    stop("signif.level must be between 0 and 1.")
  }

  # Check period lower limit
  if (!is.numeric(period.lower.limit) || period.lower.limit < 0) {
    stop("period.lower.limit must be a non-negative number.")
  }

  # --- Wavelet Transform Analysis ---

  # Store original variable and compute variance
  variable_org <- variable

  # Optional detrending before analysis
  if (detrend) {
    trend <- stats::fitted(stats::lm(variable ~ seq_along(variable)))
    variable <- variable - trend
  }

  variance1 <- stats::var(variable_org)
  n1 <- length(variable_org)

  # Standardize time series (mean=0, sd=1)
  variable <- scale(variable_org)

  # Zero-pad to next power of 2 for efficient FFT
  base2 <- floor(log2(n1) + 0.4999)
  variable <- c(variable, rep(0, (2^(base2 + 1) - n1)))
  n <- length(variable)

  # Wavelet transform parameters
  dt <- 1  # Time step (assuming unit spacing)
  dj <- 0.25  # Scale resolution (smaller = finer resolution)
  s0 <- 2 * dt  # Smallest scale
  J <- floor((1 / dj) * log2(n1 * dt / s0))  # Number of scales
  scale <- s0 * 2^((0:J) * dj)  # Scale vector

  # Wave number vector for FFT
  k <- c(0:(floor(n / 2)), -rev(1:floor((n - 1) / 2))) * ((2 * pi) / (n * dt))

  # Compute FFT of time series
  f <- stats::fft(variable)

  # Initialize wavelet transform matrix (complex)
  #wave <- array(as.complex(0), c(J + 1, n))
  wave <- matrix(complex(real = 0, imaginary = 0), nrow = J + 1, ncol = n)


  # Perform continuous wavelet transform
  params <- morlet_parameters(k0 = 6)
  fourier_factor <- params["fourier_factor"]
  coi_base <- params["coi"]
  dofmin <- params["dofmin"]


  wave <- t(sapply(seq_len(J + 1), function(a1) {
    daughter <- morlet_wavelet(k, scale[a1], k0 = 6)
    stats::fft(f * daughter, inverse = TRUE) / n
  }))


  # for (a1 in 1:(J + 1)) {
  #   # Convolve FFT with wavelet at this scale
  #   daughter <- morlet_wavelet(k, scale[a1], k0 = 6)
  #   wave[a1, ] <- stats::fft(f * daughter, inverse = TRUE) / n
  # }

  # Convert scale to period
  period <- fourier_factor * scale

  # Cone of influence (region affected by edge effects)
  #coi <- coi_base * dt * c(1e-5, 1:((n1 + 1) / 2 - 1),
  #                         rev(1:(n1 / 2 - 1)), 1e-5)

  # Fixed and clearer
  coi_indices <- c(
    1e-5,
    seq_len(floor((n1 - 1) / 2)),
    rev(seq_len(ceiling((n1 - 1) / 2))),
    1e-5
  )
  coi <- coi_base * dt * coi_indices[1:n1]

  # Trim to original length and compute power
  wave <- wave[, 1:n1, drop = FALSE]
  POWER <- abs(wave)^2

  # Global Wavelet Spectrum (time-averaged power)
  GWS <- variance1 * rowMeans(POWER)


  # --- Significance Testing ---

  # Morlet wavelet empirical parameters
  # [dofmin, Cdelta, gamma_fac, unused]
  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1]  # Degrees of freedom for Morlet
  gamma_fac <- empir[3]  # Decorrelation factor

  # Estimate lag-1 autocorrelation for noise model
  # if (noise.type == "white") {
  #   lag1 <- 0
  # } else {
  #   # Compute lag-1 autocorrelation from data
  #   lag1 <- stats::cor(variable_org[-n1], variable_org[-1])
  #   # Ensure lag1 is in valid range
  #   lag1 <- max(min(lag1, 0.999), -0.999)
  # }

  # Better approach with bias correction (Torrence & Compo 1998)
  if (noise.type == "white") {
    lag1 <- 0
  } else {
    # Bias-corrected lag-1 autocorrelation
    lag1_raw <- stats::cor(variable_org[-n1], variable_org[-1])

    # Apply bias correction for small samples
    bias_correction <- (1 + 2 * lag1_raw) / (n1 - 2)
    lag1 <- lag1_raw - bias_correction

    # Clamp to valid range
    lag1 <- max(min(lag1, 0.999), -0.999)
  }

  # Fourier power spectrum of theoretical noise
  freq <- dt / period
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2)

  # Point-wise significance for wavelet power
  chisquare <- stats::qchisq(signif.level, dofmin) / dofmin
  signif <- fft_theor * chisquare

  # Significance matrix (normalized power)
  #sigm <- POWER / (outer(signif, rep(1, n1)))
  sigm <- sweep(POWER, 1, signif, FUN = "/")


  # Effective degrees of freedom (accounts for scale decorrelation)
  dof <- n1 - scale
  dof[dof < 1] <- 1
  dof <- dofmin * sqrt(1 + (dof * dt / gamma_fac / scale)^2)
  dof[dof < dofmin] <- dofmin

  # GWS significance threshold
  chisquare_GWS <- stats::qchisq(signif.level, dof) / dof
  GWS_signif <- fft_theor * variance1 * chisquare_GWS


  # --- Identify Significant Periods ---

  # Find periods where GWS exceeds significance and period > lower limit
  sig_periods <- which(GWS > GWS_signif & period > period.lower.limit)

  # Group consecutive significant periods
  if (length(sig_periods) == 0) {
    signif_periods <- NULL
    COMPS <- NULL
  } else {
    # Split into consecutive groups
    sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))

    # For each group, select the period with maximum GWS
    signif_periods <- unlist(
      lapply(sig_periods_grp, function(x) x[which.max(GWS[x])]),
      use.names = FALSE
    )

    # Extract wavelet components for significant periods
    variable_sd <- stats::sd(variable_org)
    COMPS <- extract_wavelet_components(
      wave = wave,
      signif_periods = signif_periods,
      scale = scale,
      dj = dj,
      dt = dt,
      variable_sd = variable_sd,
      Cdelta = 0.776,
      w0_0 = pi^(-1/4)
    )

    # Add descriptive names with period values
    period_labels <- round(period[signif_periods], 2)
    names(COMPS) <- paste0("Period_", period_labels)

    # Compute residual (noise) component
    COMPS$NOISE <- variable_org - rowSums(COMPS)
  }
  # --- Return Results ---

  return(list(
    GWS = GWS,
    GWS_signif = GWS_signif,
    GWS_period = period,
    signif_periods = signif_periods,
    power = POWER,
    coi = coi,
    sigm = sigm,
    COMPS = COMPS,
    diagnostics = list(
        lag1 = lag1,                     # Estimated autocorrelation
        variance = variance1,            # Original variance
        n_original = n1,                 # Original length
        n_padded = n,                    # Padded length
        dof = dof,                       # Degrees of freedom
        scale = scale,                   # Scale vector
        fourier_factor = fourier_factor) # For user reference
  ))
}
