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

  if (!is.numeric(variable)) {
    stop("variable must be numeric")
  }

  if (anyNA(variable)) {
    stop("variable contains missing values")
  }

  if (length(variable) < 16) {
    stop("variable must have at least 16 observations")
  }

  if (!(noise.type %in% c("white", "red"))) {
    stop("noise.type must be 'white' or 'red'")
  }

  if (!is.numeric(signif.level) || length(signif.level) != 1L ||
      signif.level <= 0 || signif.level >= 1) {
    stop("signif.level must be between 0 and 1")
  }

  if (!is.numeric(period.lower.limit) || length(period.lower.limit) != 1L ||
      period.lower.limit < 0) {
    stop("period.lower.limit must be a non-negative number")
  }

  # --- Wavelet Transform Analysis ---

  variable_org <- as.numeric(variable)

  if (detrend) {
    trend <- stats::fitted(stats::lm(variable_org ~ seq_along(variable_org)))
    variable_org <- variable_org - trend
  }

  variance1 <- stats::var(variable_org)
  n1 <- length(variable_org)

  # Standardize to mean 0, sd 1 as a numeric vector (avoid scale() matrix)
  variable_std <- (variable_org - mean(variable_org)) / stats::sd(variable_org)

  # Zero-pad to next power of 2
  base2 <- floor(log2(n1) + 0.4999)
  variable_pad <- c(variable_std, rep(0, (2^(base2 + 1) - n1)))
  n <- length(variable_pad)

  dt <- 1
  dj <- 0.25
  s0 <- 2 * dt
  J <- floor((1 / dj) * log2(n1 * dt / s0))
  scale <- s0 * 2^((0:J) * dj)

  k <- c(0:(floor(n / 2)), -rev(1:floor((n - 1) / 2))) * ((2 * pi) / (n * dt))

  f <- stats::fft(variable_pad)

  params <- morlet_parameters(k0 = 6)
  fourier_factor <- params["fourier_factor"]
  coi_base <- params["coi"]

  wave <- t(sapply(seq_len(J + 1), function(a1) {
    daughter <- morlet_wavelet(k, scale[a1], k0 = 6)
    stats::fft(f * daughter, inverse = TRUE) / n
  }))

  period <- as.numeric(fourier_factor * scale)

  coi_indices <- c(
    1e-5,
    seq_len(floor((n1 - 1) / 2)),
    rev(seq_len(ceiling((n1 - 1) / 2))),
    1e-5
  )
  coi <- as.numeric(coi_base * dt * coi_indices[1:n1])

  wave <- wave[, 1:n1, drop = FALSE]
  POWER <- abs(wave)^2

  GWS <- as.numeric(variance1 * rowMeans(POWER))

  # --- Significance Testing ---

  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1]
  gamma_fac <- empir[3]

  if (noise.type == "white") {
    lag1 <- 0
  } else {
    lag1_raw <- stats::cor(variable_org[-n1], variable_org[-1])
    bias_correction <- (1 + 2 * lag1_raw) / (n1 - 2)
    lag1 <- lag1_raw - bias_correction
    lag1 <- max(min(lag1, 0.999), -0.999)
  }

  freq <- dt / period
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2)

  chisquare <- stats::qchisq(signif.level, dofmin) / dofmin
  signif <- fft_theor * chisquare

  sigm <- sweep(POWER, 1, signif, FUN = "/")

  dof <- n1 - scale
  dof[dof < 1] <- 1
  dof <- dofmin * sqrt(1 + (dof * dt / gamma_fac / scale)^2)
  dof[dof < dofmin] <- dofmin

  chisquare_GWS <- stats::qchisq(signif.level, dof) / dof
  GWS_signif <- as.numeric(fft_theor * variance1 * chisquare_GWS)

  # --- Identify Significant Periods ---

  sig_periods <- which(GWS > GWS_signif & period > period.lower.limit)

  if (length(sig_periods) == 0) {
    signif_periods <- integer(0)  # IMPORTANT: never NULL
    COMPS <- NULL
  } else {
    sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))

    signif_periods <- as.integer(unlist(
      lapply(sig_periods_grp, function(x) x[which.max(GWS[x])]),
      use.names = FALSE
    ))

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

    period_labels <- round(period[signif_periods], 2)
    names(COMPS) <- paste0("Period_", period_labels)

    COMPS$NOISE <- variable_org - rowSums(COMPS)
  }

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
      lag1 = lag1,
      variance = variance1,
      n_original = n1,
      n_padded = n,
      dof = dof,
      scale = scale,
      fourier_factor = fourier_factor
    )
  ))
}
