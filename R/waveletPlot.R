#' Plot Global Wavelet Power Spectrum Using ggplot2
#'
#' This function generates a publication-quality plot of the global wavelet power spectrum using `ggplot2`.
#' It can display observed power spectra, significance thresholds, and simulated spectra as envelopes or ribbons.
#' Designed to visually compare observed versus simulated spectral power across Fourier periods.
#'
#' @param power.period Numeric vector. Fourier periods (typically in years) returned by wavelet analysis.
#' @param power.signif Numeric vector. Significance thresholds for each period, e.g., the 95th percentile global power at each scale.
#' @param power.obs Numeric vector. Observed global wavelet power spectrum.
#' @param power.sim Optional numeric matrix. Simulated global wavelet spectra (each column = one simulation); if provided, displays simulation envelopes and mean.
#'
#' @return A `ggplot2` object representing the wavelet power spectrum plot.
#'
#' @details
#' The function provides a flexible visualization of wavelet spectra, supporting both single observed spectra and comparisons to simulated ensembles.
#' \itemize{
#'   \item When `power.sim` is `NULL`, only the observed power and significance threshold are plotted.
#'   \item When `power.sim` is provided, the plot includes:
#'     \itemize{
#'       \item The envelope (min/max) of the simulated spectra as a ribbon,
#'       \item The mean of the simulations as a solid line,
#'       \item The observed spectrum and significance threshold as lines.
#'     }
#' }
#'
#' @section Aesthetics:
#' \describe{
#'   \item{Red dashed line}{Significance threshold}
#'   \item{Blue line}{Observed power spectrum}
#'   \item{Black line}{Mean of simulated spectra (if provided)}
#'   \item{Gray ribbon}{Range of simulated spectra (min/max; if provided)}
#' }
#'
#' @examples
#' # Example data
#' period <- seq(1, 64, by = 1)
#' signif <- rep(1500, length(period))
#' obs <- 2000 * exp(-((period - 15) / 10)^2) + 1000
#' sim_mat <- sapply(1:30, function(i) obs + rnorm(length(period), sd = 200))
#'
#' # Plot observed and simulated spectra
#' p <- waveletPlot(
#'   power.period = period,
#'   power.signif = signif,
#'   power.obs = obs,
#'   power.sim = sim_mat
#' )
#' print(p)
#'
#' @import ggplot2
#' @importFrom dplyr tibble
#' @importFrom scales comma
#' @export
waveletPlot <- function(
    power.period = NULL,
    power.signif = NULL,
    power.obs = NULL,
    power.sim = NULL) {


  tsn <- length(power.period)

  if (is.null(power.sim)) {
    df <- tibble(power.period, power.signif, power.obs)

    p <- ggplot(df, aes(x = power.period)) +
      theme_light(base_size = 11) +
      geom_line(aes(y = power.signif), color = "red", linetype = "dashed", linewidth = 0.6) +
      geom_line(aes(y = power.obs), color = "blue", linewidth = 0.6)
  } else {
    savg <- apply(as.matrix(power.sim), 1, mean)
    slow <- apply(as.matrix(power.sim), 1, min)
    sup <- apply(as.matrix(power.sim), 1, max)

    # slow <- apply(power.sim, 1, function(x) quantile(x, 0.025))
    # sup  <- apply(power.sim, 1, function(x) quantile(x, 0.975))

    df <- tibble(power.period, power.signif, power.obs, slow = slow[1:tsn], sup = sup[1:tsn], savg = savg[1:tsn])

    p <- ggplot(df, aes(x = power.period)) +
      theme_light(base_size = 11) +
      geom_ribbon(aes(ymin = slow, ymax = sup), alpha = 0.2) +
      geom_line(aes(y = power.signif), color = "red", linetype = "dashed", linewidth = 0.6) +
      geom_line(aes(y = savg), color = "black", linewidth = 0.6) +
      geom_line(aes(y = power.obs), color = "blue", linewidth = 0.6)
  }

  p <- p + theme_light(base_size = 11) +
    labs(x = "Period (years)", y = expression(paste("Power (", mm^2, ")"))) #+
  # scale_x_continuous(breaks=seq(0,50,5), expand=c(0,0)) +
  # scale_y_log10(limits = c(5*10**2,5*10**5),labels = comma)

  return(p)
}
