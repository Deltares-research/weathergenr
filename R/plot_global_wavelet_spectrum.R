#' Plot Global Wavelet Spectrum
#'
#' Generates a publication quality plot of the global wavelet power spectrum.
#' The function supports comparison between an observed spectrum and an ensemble
#' of simulated spectra, displayed using an envelope and an ensemble mean.
#'
#' @param power.period Numeric vector. Fourier periods, for example in years,
#'   corresponding to the wavelet scales used in the analysis.
#' @param power.signif Numeric vector. Significance threshold for the global
#'   wavelet power at each period, for example a 95 percent confidence level.
#' @param power.obs Numeric vector. Observed global wavelet power spectrum.
#' @param power.sim Optional numeric matrix. Simulated global wavelet spectra.
#'   Rows correspond to periods and columns correspond to individual simulations.
#'   When provided, the simulation range and ensemble mean are shown.
#'
#' @return A ggplot2 object representing the global wavelet power spectrum.
#'
#' @details
#' The global wavelet spectrum summarizes wavelet power across time for each
#' Fourier period. This plot is commonly used to identify dominant variability
#' scales and to compare observed behavior against stochastic or climate driven
#' ensembles.
#'
#' Plot behavior is as follows:
#' \itemize{
#'   \item The observed global wavelet spectrum and its significance threshold
#'   are always plotted.
#'   \item When simulated spectra are provided, the ensemble range is shown
#'   as a ribbon and the ensemble mean is shown as a solid line.
#' }
#'
#' @section Aesthetics:
#' \itemize{
#'   \item Blue solid line: observed global wavelet spectrum
#'   \item Red dashed line: significance threshold
#'   \item Black solid line: ensemble mean, when provided
#'   \item Gray ribbon: ensemble envelope defined by minimum and maximum values,
#'   when provided
#' }
#'
#' @examples
#' \dontrun{
#' period <- seq(1, 64)
#' signif <- rep(1500, length(period))
#' obs <- 2000 * exp(-((period - 15) / 10)^2) + 1000
#' simmat <- sapply(
#'   1:30,
#'   function(i) obs + rnorm(length(period), sd = 200)
#' )
#'
#' p <- plot_global_wavelet_spectrum(
#'   power.period = period,
#'   power.signif = signif,
#'   power.obs = obs,
#'   power.sim = simmat
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

  # if (length(power.signif) != n_period ||
  #     length(power.obs)    != n_period) {
  #   stop(
  #     "power.period, power.signif, and power.obs must have identical length.",
  #     call. = FALSE
  #   )
  # }

  if (!is.null(power.sim)) {
    power.sim <- as.matrix(power.sim)
  }

  # -------------------------------------------------------------------------
  # Construct plotting data
  # -------------------------------------------------------------------------



  obs_length <- length(power.obs)
  sim_length <- nrow(power.sim)

  df_obs <- data.frame(period = power.period[1:obs_length], obs = power.obs,
    signif = power.signif[1:obs_length], stringsAsFactors = FALSE)

  # Base plot with observed spectrum and significance threshold
  p <- ggplot(df_obs, aes(x = period)) +
    theme_light(base_size = 11) +
    geom_line(aes(y = obs), color = "blue", linewidth = 0.6) +
    geom_line(aes(y = signif), color = "red", linetype = "dashed",
              linewidth = 0.6, na.rm = TRUE)


  # Add simulation envelope if provided
  if (!is.null(power.sim)) {


    sim_mean <- rowMeans(power.sim, na.rm = TRUE)
    sim_min  <- apply(power.sim, 1, min, na.rm = TRUE)
    sim_max  <- apply(power.sim, 1, max, na.rm = TRUE)

    # Create data frame for simulation data
    df_sim <- tibble(period = power.period, signif = power.signif,
        sim = power.sim, sim_lo = sim_min, sim_hi = sim_max,
        sim_mu = sim_mean)

    # Add simulation ribbon and mean line
    p <- p +
      #geom_ribbon(aes(x = period, ymin = sim_lo, ymax = sim_hi), data = df_sim,
      #            fill = "grey60", alpha = 0.25, na.rm = TRUE) +
      geom_line(aes(x = period, y = sim_mu), data = df_sim,
                color = "black", linewidth = 0.6, na.rm = TRUE)
  }

  # Add scales and labels after all data layers
  p <- p +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Period (years)", y = expression(paste("Power (", mm^2, ")")))

  p
}
