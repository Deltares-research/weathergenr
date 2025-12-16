#' Plot Global Wavelet Spectrum
#'
#' Generate a publication-quality plot of the global wavelet power spectrum (GWS).
#' The function supports comparison between an observed spectrum and an ensemble
#' of simulated spectra, displayed as an envelope and mean.
#'
#' @param power.period Numeric vector. Fourier periods (e.g., years) corresponding
#'   to the wavelet scales used in the analysis.
#' @param power.signif Numeric vector. Significance threshold for the global wavelet
#'   power at each period (e.g., 95% confidence level).
#' @param power.obs Numeric vector. Observed global wavelet power spectrum.
#' @param power.sim Optional numeric matrix. Simulated global wavelet spectra with
#'   dimensions `[n_periods, n_simulations]`. When provided, the simulation range
#'   (min-max) and ensemble mean are shown.
#'
#' @return A `ggplot2` object representing the global wavelet power spectrum.
#'
#' @details
#' The global wavelet spectrum (GWS) summarizes wavelet power across time for each
#' Fourier period. This plot is typically used to identify dominant variability
#' scales and to compare observed behavior against stochastic or climate-driven
#' ensembles.
#'
#' Plot behavior:
#'
#' - Always plots the observed GWS and its significance threshold.
#' - When `power.sim` is provided:
#'   - The ensemble envelope (min-max) is shown as a ribbon.
#'   - The ensemble mean is shown as a solid line.
#'
#' @section Aesthetics:
#'
#' - **Blue solid line**: Observed global wavelet spectrum
#' - **Red dashed line**: Significance threshold
#' - **Black solid line**: Ensemble mean (if provided)
#' - **Gray ribbon**: Ensemble envelope (min-max; if provided)
#'
#' @examples
#' \dontrun{
#' period <- seq(1, 64)
#' signif <- rep(1500, length(period))
#' obs    <- 2000 * exp(-((period - 15) / 10)^2) + 1000
#' simmat <- sapply(1:30, function(i) obs + rnorm(length(period), sd = 200))
#'
#' p <- plot_global_wavelet_spectrum(
#'   power.period = period,
#'   power.signif = signif,
#'   power.obs    = obs,
#'   power.sim    = simmat
#' )
#' print(p)
#' }
#'
#' @import ggplot2
#' @importFrom tibble tibble
#' @export
plot_global_wavelet_spectrum <- function(
    power.period,
    power.signif,
    power.obs,
    power.sim = NULL
) {

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  stopifnot(
    is.numeric(power.period),
    is.numeric(power.signif),
    is.numeric(power.obs)
  )

  n_period <- length(power.period)

  if (length(power.signif) != n_period ||
      length(power.obs)    != n_period) {
    stop(
      "power.period, power.signif, and power.obs must have identical length.",
      call. = FALSE
    )
  }

  if (!is.null(power.sim)) {
    power.sim <- as.matrix(power.sim)
    if (nrow(power.sim) != n_period) {
      stop(
        "Rows of power.sim must equal length of power.period.",
        call. = FALSE
      )
    }
  }

  # -------------------------------------------------------------------------
  # Construct plotting data
  # -------------------------------------------------------------------------
  if (is.null(power.sim)) {

    df <- tibble(
      period = power.period,
      signif = power.signif,
      obs    = power.obs
    )

    p <- ggplot(df, aes(x = period)) +
      geom_line(aes(y = signif),
                color = "red", linetype = "dashed", linewidth = 0.6) +
      geom_line(aes(y = obs),
                color = "blue", linewidth = 0.6)

  } else {

    sim_mean <- rowMeans(power.sim, na.rm = TRUE)
    sim_min  <- apply(power.sim, 1, min, na.rm = TRUE)
    sim_max  <- apply(power.sim, 1, max, na.rm = TRUE)

    df <- tibble(
      period = power.period,
      signif = power.signif,
      obs    = power.obs,
      sim_lo = sim_min,
      sim_hi = sim_max,
      sim_mu = sim_mean
    )

    p <- ggplot(df, aes(x = period)) +
      geom_ribbon(aes(ymin = sim_lo, ymax = sim_hi),
                  fill = "grey60", alpha = 0.25) +
      geom_line(aes(y = signif),
                color = "red", linetype = "dashed", linewidth = 0.6) +
      geom_line(aes(y = sim_mu),
                color = "black", linewidth = 0.6) +
      geom_line(aes(y = obs),
                color = "blue", linewidth = 0.6)
  }

  # -------------------------------------------------------------------------
  # Common styling
  # -------------------------------------------------------------------------
  p +
    theme_light(base_size = 11) +
    labs(
      x = "Period (years)",
      y = expression(paste("Power (", mm^2, ")"))
    )
}
