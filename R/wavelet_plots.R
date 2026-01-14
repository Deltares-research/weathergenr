#' Plot Wavelet Power and Global Wavelet Spectrum
#'
#' @description
#' Produces a two-panel diagnostic figure for a 1D time series wavelet analysis:
#' (1) the time-period wavelet power field and (2) the global wavelet spectrum
#' (time-averaged power by period) with its significance threshold.
#'
#' @details
#' \strong{What this plot is for}
#' \itemize{
#'   \item Identify when specific periodicities are active in time (power field).
#'   \item Summarize dominant variability scales across the full record (global spectrum).
#'   \item Check edge effects and statistical significance (COI and contour).
#' }
#'
#' \strong{Key plotting choices}
#' \itemize{
#'   \item Power is shown on a \eqn{\log_2} scale after discarding non-finite values and
#'   non-positive power.
#'   \item The power field is drawn on the native wavelet grid (no interpolation) using
#'   explicit cell boundaries, which avoids smoothing artifacts.
#'   \item Color limits are set using robust quantiles (5th-95th percentile) to avoid
#'   low-contrast fields when dynamic range is large.
#' }
#'
#' \strong{Inputs must be consistent}
#' \itemize{
#'   \item \code{n_time = length(series)} and \code{n_period = length(period)}.
#'   \item \code{power} and \code{signif_mask} must have dimensions
#'   \code{n_period x n_time}.
#'   \item \code{coi} must have length \code{n_time}.
#'   \item \code{period} must be strictly increasing.
#'   \item If supplied, \code{time} must be strictly increasing and have length \code{n_time}.
#' }
#'
#' @param series Numeric vector. The time series used to define the time dimension
#'   (length must match the number of columns in \code{power}).
#' @param time Optional numeric vector of length \code{length(series)} used for the x-axis.
#'   Must be strictly increasing. If \code{NULL}, uses \code{seq_len(length(series))}.
#' @param period Numeric vector of wavelet periods (scales), typically in years (or days),
#'   strictly increasing.
#' @param power Numeric matrix of wavelet power with dimensions
#'   \code{length(period) x length(series)}.
#' @param gws Numeric vector. Global wavelet spectrum (mean power for each \code{period}).
#'   Must have length \code{length(period)}.
#' @param gws_signif Numeric vector. Significance threshold for \code{gws}. Must have length
#'   \code{length(period)}.
#' @param coi Numeric vector of length \code{length(series)} giving the cone of influence
#'   in the same units as \code{period}. Values outside the plotted period range are ignored.
#' @param signif_mask Numeric matrix with dimensions \code{length(period) x length(series)}.
#'   A contour is drawn at \code{1}. This is typically a significance ratio or a binary mask
#'   encoded as 0/1 (or NA).
#' @param unit Character. Unit label used in the global spectrum axis label (default: \code{"mm"}).
#'
#' @return A \code{patchwork} object combining the power-field panel and the global spectrum panel.
#'
#' @seealso
#' \code{\link{plot_wavelet_global_spectrum}} for a standalone global spectrum plot.
#'
#' @examples
#' \dontrun{
#' # Assuming you computed wavelet outputs elsewhere:
#' # w <- analyze_wavelet_spectrum(series)
#' p <- plot_wavelet_power(
#'   series      = series,
#'   time        = years,
#'   period      = w$period,
#'   power       = w$power,
#'   gws         = w$gws,
#'   gws_signif  = w$gws_signif,
#'   coi         = w$coi,
#'   signif_mask = w$signif_mask,
#'   unit        = "mm"
#' )
#' print(p)
#' }
#'
#' @export
#' @import ggplot2
#' @import patchwork
#' @import scales
#' @importFrom stats quantile
plot_wavelet_power <- function(
    series,
    time = NULL,
    period,
    power,
    gws,
    gws_signif,
    coi,
    signif_mask,
    unit = "mm"
) {

  if (!is.numeric(series)) stop("'series' must be numeric.", call. = FALSE)
  n_time <- length(series)
  n_period <- length(period)

  if (!is.matrix(power) || nrow(power) != n_period || ncol(power) != n_time) {
    stop("Invalid 'power' dimensions.", call. = FALSE)
  }
  if (!is.matrix(signif_mask)) stop("'signif_mask' must be a matrix.", call. = FALSE)

  time_axis <- if (is.null(time)) seq_len(n_time) else time

  .cell_bounds <- function(x) {
    mids <- (x[-1L] + x[-length(x)]) / 2
    list(
      lower = c(x[1L] - (mids[1L] - x[1L]), mids),
      upper = c(mids, x[length(x)] + (x[length(x)] - mids[length(mids)]))
    )
  }

  xb <- .cell_bounds(time_axis)
  yb <- .cell_bounds(period)

  power[!is.finite(power) | power <= 0] <- NA_real_
  z <- log2(power)

  df_power <- data.frame(
    x = rep(time_axis, times = n_period),
    y = rep(period, each = n_time),
    xmin = rep(xb$lower, times = n_period),
    xmax = rep(xb$upper, times = n_period),
    ymin = rep(yb$lower, each = n_time),
    ymax = rep(yb$upper, each = n_time),
    z = as.vector(t(z))
  )
  df_power <- df_power[is.finite(df_power$z), ]

  zlims <- quantile(df_power$z, c(0.05, 0.95), na.rm = TRUE)

  p_spectrum <- ggplot(df_power) +
    theme_light() +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = z)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(limits = zlims, oob = squish) +
    labs(x = "Time", y = "Period (years)") +
    guides(fill = "none") +
    geom_line(data = data.frame(x = time_axis, y = coi),
              aes(x = x, y = y),
              linetype = "dashed", color = "red") +
    stat_contour(
      data = data.frame(
        x = rep(time_axis, times = n_period),
        y = rep(period, each = n_time),
        z = as.vector(t(signif_mask))
      ),
      aes(x = x, y = y, z = z),
      breaks = 1,
      color = "black"
    )

  gws_df <- data.frame(period = period, gws = gws, signif = gws_signif)

  p_gws <- ggplot(gws_df, aes(x = period, y = gws)) +
    theme_light() +
    geom_line() +
    geom_point(size = 1) +
    geom_line(aes(y = signif), color = "red", linetype = "dashed") +
    coord_flip() +
    scale_x_reverse(limits = range(period)) +
    labs(y = bquote(Power ~ (.(unit)^2)), x = "")

  p_spectrum + p_gws + plot_layout(widths = c(3, 1.2))
}


#' Plot Global Wavelet Spectrum
#'
#' @description
#' Plots the global wavelet spectrum (time-averaged wavelet power by period)
#' for an observed series, together with a period-wise significance threshold.
#' Optionally overlays the ensemble mean of simulated global spectra.
#'
#' @details
#' \strong{What this plot is for}
#' \itemize{
#'   \item Identify dominant variability scales (peaks in global power).
#'   \item Compare observed spectral structure to a simulated ensemble.
#'   \item Interpret variability relative to a significance threshold.
#' }
#'
#' \strong{Ensemble behavior}
#' When \code{sim_power} is provided, the function adds the ensemble mean as a solid line.
#' This function intentionally keeps the default view simple; if you want an ensemble
#' envelope (e.g., min/max or quantiles) add it outside this function.
#'
#' @param period Numeric vector of wavelet periods (scales), typically in years (or days).
#'   Should be strictly increasing for meaningful interpretation.
#' @param signif Numeric vector of length \code{length(period)} giving the global-spectrum
#'   significance threshold at each period.
#' @param obs_power Numeric vector of length \code{length(period)} giving the observed
#'   global wavelet spectrum.
#' @param sim_power Optional numeric matrix of simulated global spectra with
#'   \code{nrow(sim_power) == length(period)} and one column per simulation.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso
#' \code{\link{plot_wavelet_power}} for a combined diagnostic plot including the time-period power field.
#'
#' @examples
#' \dontrun{
#' p <- plot_wavelet_global_spectrum(
#'   period   = w$period,
#'   signif   = w$gws_signif,
#'   obs_power = w$gws,
#'   sim_power = w$gws_sim  # matrix: period x n_sim (optional)
#' )
#' print(p)
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom tibble tibble
plot_wavelet_global_spectrum <- function(
    period,
    signif,
    obs_power,
    sim_power = NULL
) {

  df_obs <- data.frame(
    period = period,
    obs = obs_power,
    signif = signif
  )

  p <- ggplot(df_obs, aes(x = period)) +
    theme_light() +
    geom_line(aes(y = obs), color = "blue") +
    geom_line(aes(y = signif), color = "red", linetype = "dashed") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Period (years)", y = expression(paste("Power (", mm^2, ")")))

  if (!is.null(sim_power)) {
    sim_mu <- rowMeans(sim_power, na.rm = TRUE)
    p <- p + geom_line(
      data = tibble(period = period, sim_mu = sim_mu),
      aes(x = period, y = sim_mu),
      color = "black"
    )
  }

  p
}
