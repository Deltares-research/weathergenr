#' Wavelet Decomposition of Time Series into Low-Frequency Components and Noise
#'
#' @description
#' Decomposes a univariate time series into significant low-frequency (periodic) components and residual noise
#' using Morlet wavelet analysis, based on periods identified as significant by prior wavelet analysis.
#' Each significant periodic component is extracted and returned as a column, along with the reconstructed "noise" (high-frequency residual).
#'
#' @details
#' This function is used in the Wavelet Autoregressive Modeling (WARM) framework to isolate low-frequency signals
#' (such as decadal or multi-annual climate cycles) for separate stochastic modeling, while treating the remaining noise
#' as a stationary process. The decomposition uses a Morlet wavelet transform and reconstructs time series components
#' at each significant period (as found by prior `waveletAnalysis`), plus a residual ("noise") component.
#'
#' @param variable Numeric vector. The time series to be decomposed (e.g., annual precipitation).
#' @param signif.periods List. Each element is a vector of integer scale indices representing significant periods,
#'    as returned by \code{waveletAnalysis()} (see its \code{signif_periods} output).
#' @param noise.type Character string. Background noise type: either `"white"` (default) or `"red"`.
#' @param signif.level Numeric between 0 and 1. Significance level for wavelet significance testing (default 0.90).
#' @param plot Logical. If TRUE, saves a plot of decomposition (default TRUE).
#' @param output.path String. Directory to save the plot if `plot=TRUE` (default: NULL).
#'
#' @import dplyr tidyr ggplot2
#' @return
#' A tibble (data frame) with columns:
#'   \item{Component_1, Component_2, ...}{Numeric vectors: extracted low-frequency signals for each significant period.}
#'   \item{NOISE}{Numeric vector: residual high-frequency noise.}
#'
#' @seealso
#'   \code{\link{waveletAnalysis}} for identification of significant periods.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' x <- as.numeric(scale(rnorm(50))) # Example annual series
#' sig_periods <- list(c(3, 5)) # Assume periods 3 and 5 are significant (from waveletAnalysis)
#' comp <- waveletDecompose(
#'   variable = x,
#'   signif.periods = sig_periods,
#'   plot = FALSE
#' )
#' head(comp)
#' }
#'
#' @export
waveletDecompose <- function(variable = NULL,
                             signif.periods = NULL,
                             noise.type = "white",
                             signif.level = 0.90,
                             plot = TRUE,
                             output.path = NULL) {
  if (is.null(signif.periods)) {
    stop("signif.periods must not be NULL")
  }
  if (length(signif.periods) == 0) {
    # If empty, return just the noise (no components)
    return(tibble::tibble(NOISE = variable))
  }

  # Workaround for rlang warning
  value <- year <- 0

  # Number of orthogonal component series that representing a low-freq signal
  NUM_FINAL_PERIODS <- length(signif.periods) # OK

  # List of all significant periods c(2,3,4...)
  ALL_SIG_PERIODS <- unlist(signif.periods, use.names = FALSE) # OK

  # Length of all components (this needs to be length of components e.g. c(len(COMP1), leng(COMP2)
  NUM_PERIODS_ALL_COMPS <- as.numeric(sapply(signif.periods, function(x) length(x))) # OK

  #### Define Wavelet function used that returns transform
  waveletf <- function(k, s) {
    nn <- length(k)
    k0 <- 6 # nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
    z <- array(1, nn)
    z[which(k <= 0)] <- 0
    expnt <- -((s * k - k0)^2 / 2) * z
    norm <- sqrt(s * k[2]) * (pi^(-0.25)) * sqrt(nn) # total energy=N   [Eqn(7)]
    daughter <- norm * exp(expnt)
    daughter <- daughter * z
    fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
    coi <- fourier_factor / sqrt(2) # Cone-of-influence [Sec.3g]
    dofmin <- 2 # Degrees of freedom
    return(daughter)
  }

  # Define Wavelet function used that returns fourier_factor,
  waveletf2 <- function(k, s) {
    nn <- length(k)
    k0 <- 6 # nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
    z <- array(1, nn)
    z[which(k <= 0)] <- 0
    expnt <- -((s * k - k0)^2 / 2) * z
    norm <- sqrt(s * k[2]) * (pi^(-0.25)) * sqrt(nn) # total energy=N   [Eqn(7)]
    daughter <- norm * exp(expnt)
    daughter <- daughter * z
    fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
    coi <- fourier_factor / sqrt(2) # Cone-of-influence [Sec.3g]
    dofmin <- 2 # Degrees of freedom
    return(c(fourier_factor, coi, dofmin))
  }

  # :::::::::::::::: PERFORM WAVELET   ::::::::::::::::::::::::::::::::::::::::::

  # ....construct time series to analyze, pad if necessary
  current_variable_org <- variable
  variance1 <- stats::var(current_variable_org)
  n1 <- length(current_variable_org)

  current_variable <- scale(current_variable_org)
  variance2 <- stats::var(current_variable)
  base2 <- floor(log(n1) / log(2) + 0.4999) # power of 2 nearest to N

  current_variable <- c(current_variable, rep(0, (2^(base2 + 1) - n1)))
  n <- length(current_variable)

  # Determine parameters for Wavelet analysis
  dt <- 1
  dj <- 0.25
  s0 <- 2 * dt
  J <- floor((1 / dj) * log((n1 * dt / s0), base = 2))

  # ....construct SCALE array & empty PERIOD & WAVE arrays
  scale <- s0 * 2^((0:J) * dj)
  period <- scale
  wave <- array(as.complex(0), c(J + 1, n)) # define the wavelet array

  # ....construct wavenumber array used in transform [Eqn(5)]
  k <- c(1:floor(n / 2))
  k <- k * ((2. * pi) / (n * dt))
  k <- c(0, k, -rev(k[1:floor((n - 1) / 2)]))

  # fourier transform of standardized precipitation
  f <- stats::fft(current_variable, inverse = FALSE)

  # loop through all scales and compute transform
  for (a1 in 1:(J + 1)) {
    daughter <- waveletf(k, scale[a1])
    results <- waveletf2(k, scale[a1])
    fourier_factor <- results[1]
    coi <- results[2]
    dofmin <- results[3]
    wave[a1, ] <- stats::fft(f * daughter, inverse = TRUE) / n # wavelet transform[Eqn(4)]
  }

  period <- fourier_factor * scale
  coi <- coi * dt * c((0.00001), 1:((n1 + 1) / 2 - 1), rev((1:(n1 / 2 - 1))), (0.00001)) # COI [Sec.3g]
  wave <- wave[, 1:n1] # get rid of padding before returning
  POWER <- abs(wave)^2
  GWS <- variance1 * apply(POWER, FUN = mean, c(1)) # Global Wavelet Spectrum

  # :::::::::::::::: Signficance Testing  :::::::::::::::::::::::::::::::::::::::

  # get the appropriate parameters [see Table(2)]
  k0 <- 6
  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1] # Degrees of freedom with no smoothing
  Cdelta <- empir[2] # reconstruction factor
  gamma_fac <- empir[3] # time-decorrelation factor
  dj0 <- empir[4] # scale-decorrelation factor

  # for red noise background, lag1 autocorrelation = 0.72,
  # for white noise background, lag1 autocorrelation = 0
  if (noise.type == "white") {
    lag1 <- 0
  } else {
    lag1 <- .72
  }

  freq <- dt / period # normalized frequency
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2) # [Eqn(16)]
  fft_theor <- fft_theor # include time-series variance
  dof <- dofmin

  # ENTIRE POWER SPECTRUM
  chisquare <- stats::qchisq(signif.level, dof) / dof
  signif <- fft_theor * chisquare # [Eqn(18)]
  sig95 <- ((signif)) %o% (array(1, n1)) # expand signif --> (J+1)x(N) array
  sig95 <- POWER / sig95 # where ratio > 1, power is significant

  # TIME_AVERAGED (GLOBAL WAVELET SPECTRUM)
  dof <- n1 - scale
  if (length(dof) == 1) {
    dof <- array(0, (J + 1)) + dof
  }
  dof[which(dof < 1)] <- 1
  dof <- dofmin * sqrt(1 + (dof * dt / gamma_fac / scale)^2)
  tt <- which(dof < dofmin)
  dof[tt] <- dofmin
  chisquare_GWS <- array(NA, (J + 1))
  signif_GWS <- array(NA, (J + 1))
  for (a1 in 1:(J + 1)) {
    chisquare_GWS[a1] <- stats::qchisq(signif.level, dof[a1]) / dof[a1]
    signif_GWS[a1] <- fft_theor[a1] * variance1 * chisquare_GWS[a1]
  }

  # :::::::::::::::: DEFINE TIME-SERIES COMPONENTS  :::::::::::::::::::::::::::::

  COMPS <- array(0, c(length(variable), NUM_FINAL_PERIODS)) %>%
    tibble::as_tibble(.name_repair = ~ paste0("COMP", 1:NUM_FINAL_PERIODS))

  # Loop through each component
  for (i in 1:NUM_FINAL_PERIODS) {
    CUR_PERIODS <- ALL_SIG_PERIODS[1:NUM_PERIODS_ALL_COMPS[i]]
    if (i > 1) {
      CUR_PERIODS <- ALL_SIG_PERIODS[(1 + (i - 1) * NUM_PERIODS_ALL_COMPS[i - 1]):(NUM_PERIODS_ALL_COMPS[i] + (i - 1) * NUM_PERIODS_ALL_COMPS[i - 1])]
    }

    sj <- scale[CUR_PERIODS]

    # for Morlet Wavelet with freq = 6
    Cdelta <- .776
    w0_0 <- pi^(-1 / 4)
    if (length(CUR_PERIODS) > 1) {
      COMPS[, i] <- apply(stats::sd(current_variable_org) * (dj * sqrt(dt) / (Cdelta * w0_0)) * Re(wave)[CUR_PERIODS, ] / sqrt(sj), FUN = sum, c(2))
    }
    if (length(CUR_PERIODS) == 1) {
      COMPS[, i] <- stats::sd(current_variable_org) * (dj * sqrt(dt) / (Cdelta * w0_0)) * Re(wave)[CUR_PERIODS, ] / sqrt(sj)
    }
  }

  names(COMPS) <- paste0("Component_", 1:NUM_FINAL_PERIODS)
  Noise <- current_variable_org - apply(COMPS, 1, sum)

  if (plot == TRUE) {
    df <- tibble(year = 1:length(variable), Original = variable) %>%
      bind_cols(COMPS) %>%
      bind_cols(Noise = Noise) %>%
      tidyr::gather(key = variable, value = value, -year) %>%
      dplyr::mutate(variable = factor(variable,
        levels = c("Original", names(COMPS), "Noise")
      ))


    p <- ggplot2::ggplot(df, aes(x = year, y = value)) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::facet_wrap(~variable, ncol = 1, scales = "free") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Time (year)", y = "")

    # Save plot to file
    ggplot2::ggsave(file.path(output.path, "warm_decomposition.png"), height = 6, width = 8)
  }

  return(COMPS %>% dplyr::mutate(NOISE = Noise))
}
