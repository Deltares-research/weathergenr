#' Wavelet Spectral Analysis with Significance Testing
#'
#' Performs continuous wavelet transform analysis on a time series,
#' computes the global wavelet spectrum (GWS) and significance levels,
#' and (optionally) produces a plot of power spectrum.
#'
#' @param variable Numeric vector, the time series (e.g. annual precipitation).
#' @param signif.level Numeric between 0 and 1. Significance level for test (default = 0.90).
#' @param noise.type "white" (default) or "red". Background noise for spectrum test.
#' @param variable.unit Character; unit for variable (for plotting).
#' @param plot Logical, whether to save a plot (default FALSE).
#' @param output.path Optional string. Directory where plot will be saved (if plot = TRUE).
#'
#' @return A list with:
#'   \item{GWS}{Global Wavelet Spectrum}
#'   \item{GWS_signif}{Significance threshold for GWS}
#'   \item{GWS_period}{Period axis for GWS}
#'   \item{signif_periods}{Indices of significant periods}
#'
#' @examples
#' set.seed(123)
#' # Simulate AR(1) process with periodicity at 8 years
#' years <- 1:64
#' signal <- sin(2 * pi * years / 8) + arima.sim(n = 64, model = list(ar = 0.7))
#' res <- waveletAnalysis(signal, plot = FALSE)
#' print(res$GWS_period[res$signif_periods])
#'
#' @import patchwork
#' @import dplyr
#' @export
waveletAnalysis <- function(variable,
                            signif.level = 0.90,
                            noise.type = "white",
                            variable.unit = "mm",
                            plot = FALSE,
                            output.path = NULL
) {


  # --- Input Valiation
  stopifnot(is.numeric(variable), length(variable) > 8)
  if (anyNA(variable)) stop("Variable contains missing values.")
  if (!(noise.type %in% c("white", "red"))) stop("noise.type must be 'white' or 'red'")
  if (!is.numeric(signif.level) || signif.level <= 0 || signif.level >= 1) {
    stop("signif.level must be between 0 and 1.")
  }

  ##############################################################################

  #### Wavelet decomposition helper function
  extract_wavelet_components <- function(wave, signif.periods, scale, dj = 0.25, dt = 1,
        Cdelta = 0.776,  w0_0 = pi^(-1 / 4), variable_sd) {

    num_periods <- length(signif.periods)
    n <- ncol(wave)
    COMPS <- matrix(0, nrow = n, ncol = num_periods)

    for (i in seq_len(num_periods)) {

      cur_periods <- signif.periods[[i]]
      sj <- scale[cur_periods]

      # Extract real part of the wavelet coefficients for the selected periods
      W <- Re(wave)[cur_periods, , drop = FALSE]
      # Reconstruction factor (constant for this variable)
      recon_fac <- variable_sd * (dj * sqrt(dt) / (Cdelta * w0_0))
      # Divide each row by sqrt(sj)
      W_scaled <- sweep(W, 1, sqrt(sj), "/")
      # Sum across selected periods if more than one
      if (nrow(W_scaled) > 1) {
        component <- recon_fac * colSums(W_scaled)
      } else {
        component <- as.numeric(recon_fac * W_scaled)
      }
      COMPS[, i] <- component
    }
    # Convert to tibble and name columns
    comps_tb <- tibble::as_tibble(COMPS, .name_repair = ~paste0("Component_", seq_len(num_periods)))
    return(comps_tb)
  }

  # Helper: Morlet wavelet (Fourier domain)
  waveletf <- function(k, s) {
    nn <- length(k)
    k0 <- 6
    z <- as.numeric(k > 0)
    expnt <- -((s * k - k0)^2 / 2) * z
    norm <- sqrt(s * k[2]) * (pi^(-0.25)) * sqrt(nn)
    daughter <- norm * exp(expnt) * z
    return(daughter)
  }

  # Helper: Morlet parameters (for period/coi/dofmin)
  waveletf2 <- function(k, s) {
    k0 <- 6
    fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0^2))
    coi <- fourier_factor / sqrt(2)
    dofmin <- 2
    c(fourier_factor, coi, dofmin)
  }

  ##############################################################################


  # ---- Wavelet Analysis ----
  # Perform wavelet decomposition and calculate power spectrum)

  # # Standardize and zero-pad the time series
  variable_org <- variable
  variance1 <- stats::var(variable_org)
  n1 <- length(variable_org)
  variable <- scale(variable_org)
  base2 <- floor(log2(n1) + 0.4999)
  variable <- c(variable, rep(0, (2^(base2 + 1) - n1)))
  n <- length(variable)

  # Wavelet transform parameters
  dt <- 1
  dj <- 0.25
  s0 <- 2 * dt
  J <- floor((1 / dj) * log((n1 * dt / s0), base = 2))
  scale <- s0 * 2^((0:J) * dj)
  k <- c(0:(floor(n / 2)), -rev(1:floor((n - 1) / 2))) * ((2 * pi) / (n * dt))


  # Compute FFT and wavelet transform
  f <- stats::fft(variable)
  wave <- array(as.complex(0), c(J + 1, n))
  for (a1 in 1:(J + 1)) {
    daughter <- waveletf(k, scale[a1])
    wave[a1, ] <- stats::fft(f * daughter, inverse = TRUE) / n
    if (a1 == 1) {
      params <- waveletf2(k, scale[a1])
      fourier_factor <- params[1]
      coi_base <- params[2]
    }
  }
  period <- fourier_factor * scale
  coi <- coi_base * dt * c(0.00001, 1:((n1 + 1) / 2 - 1), rev(1:(n1 / 2 - 1)), 0.00001)
  wave <- wave[, 1:n1]
  POWER <- abs(wave)^2
  GWS <- variance1 * rowMeans(POWER)

  # ---- Significance Testing ----
  # Assess significance of observed power vs. noise background)

  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1]
  gamma_fac <- empir[3]
  lag1 <- if (noise.type == "white") 0 else 0.72
  freq <- dt / period
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2)
  chisquare <- stats::qchisq(signif.level, dofmin) / dofmin
  signif <- fft_theor * chisquare
  sigm <- POWER / (outer(signif, rep(1, n1)))
  dof <- n1 - scale
  dof[dof < 1] <- 1
  dof <- dofmin * sqrt(1 + (dof * dt / gamma_fac / scale)^2)
  dof[dof < dofmin] <- dofmin
  chisquare_GWS <- stats::qchisq(signif.level, dof) / dof
  GWS_signif <- fft_theor * variance1 * chisquare_GWS

  # Identify periods where GWS exceeds significance threshold
  period_lower_limit <- 0
  sig_periods <- which(GWS > GWS_signif & period > period_lower_limit)
  sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))
  signif_periods <- unlist(lapply(sig_periods_grp,
                                  function(x) as.integer(round(stats::median(x)))), use.names = FALSE)

  # If no low-frequency signal present, exit and return WARM outputs
  if(any(is.na(signif_periods))) {

    return(list(
      GWS = GWS,
      GWS_signif = GWS_signif,
      GWS_period = period,
      signif_periods = signif_periods))

  } else {

    # EXTRACT WAVELET COMPONENTS -------------------------------------------------
    variable_sd <- stats::sd(variable_org)

    COMPS <- extract_wavelet_components(
      wave = wave,
      signif.periods = signif_periods,
      scale = scale,
      dj = dj,
      dt = dt,
      variable_sd = variable_sd
    )

    names(COMPS) <- paste0("Component_", 1:length(signif_periods))
    COMPS$NOISE = variable_org - apply(COMPS, 1, sum)

    ##############################################################################

    # --- Plotting (optional)
    if (plot) {

      if (is.null(output.path)) output.path <- tempdir()
      if (!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)

      p <- plot_wavelet_spectra(
        variable = variable_org,
        period = period,
        POWER = POWER,
        GWS = GWS,
        GWS_signif = GWS_signif,
        coi = coi,
        sigm = sigm)

      ggsave(file.path(output.path, "warm_hist_analysis.png"),
             width = 8, height = 6)
    }

    # --- Remove signif periods that are > 1/2 of the time-series length


    signif_periods_final <- signif_periods[which(period[signif_periods] < floor(n1/2))]

    return(list(
      GWS = GWS,
      GWS_signif = GWS_signif,
      GWS_period = period,
      signif_periods = signif_periods_final,
      COMPS = COMPS))

  }


}

