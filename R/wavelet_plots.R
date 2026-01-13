#' Plot Wavelet Power Spectrum and Global Wavelet Spectrum
#'
#' Creates a two-panel figure showing (i) the time-period wavelet power spectrum
#' (with cone of influence and significance contour) and (ii) the global wavelet
#' spectrum for a time series.
#'
#' Key plotting choices for interpretability:
#' \itemize{
#'   \item Power is plotted on a log2 scale with robust color limits (5th-95th percentiles)
#'         to avoid a "washed out" field when the dynamic range is large.
#'   \item The spectrum is plotted on the native wavelet grid (no interpolation), avoiding
#'         smoothing/extrapolation artifacts that can flatten contrast.
#'   \item Uses a perceptually uniform \code{viridis} color scale.
#' }
#'
#' @param variable Numeric vector. The original time series (for length).
#' @param variable.year Numeric vector. Year labels for x-axis (optional).
#' @param period Numeric vector. Periods (scales) used in the wavelet transform.
#' @param POWER Matrix. Wavelet power spectrum (periods x time).
#' @param GWS Numeric vector. Global wavelet spectrum (mean power at each period).
#' @param GWS_signif Numeric vector. Significance threshold for the global wavelet spectrum.
#' @param coi Numeric vector. Cone of influence for each time point.
#' @param sigm Matrix. Significance mask or ratio (periods x time). Contour is drawn at 1.
#' @param variable.unit Character. Unit label for the variable (default is "mm").
#'
#' @importFrom stats var sd cor fft simulate
#' @return A patchwork plot object with the time-period power spectrum and global wavelet spectrum.
#' @export
plot_wavelet_spectra <- function(
    variable,
    variable.year = NULL,
    period,
    POWER,
    GWS,
    GWS_signif,
    coi,
    sigm,
    variable.unit = "mm") {

  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (is.null(variable) || is.null(period) || is.null(POWER) ||
      is.null(GWS) || is.null(GWS_signif) || is.null(coi) || is.null(sigm)) {
    stop("All input arguments must be provided and non-null.", call. = FALSE)
  }
  if (!is.numeric(variable) || !is.vector(variable)) {
    stop("'variable' must be a numeric vector.", call. = FALSE)
  }

  n_time <- length(variable)

  if (!is.numeric(period) || !is.vector(period) || length(period) < 2L) {
    stop("'period' must be a numeric vector of length >= 2.", call. = FALSE)
  }
  n_period <- length(period)

  if (!is.matrix(POWER) || nrow(POWER) != n_period || ncol(POWER) != n_time) {
    stop("POWER must be a matrix with dim = length(period) x length(variable).", call. = FALSE)
  }
  if (!is.matrix(sigm) || nrow(sigm) != n_period || ncol(sigm) != n_time) {
    stop("sigm must be a matrix with dim = length(period) x length(variable).", call. = FALSE)
  }
  if (!is.numeric(GWS) || length(GWS) != n_period) {
    stop("GWS must have length equal to length(period).", call. = FALSE)
  }
  if (!is.numeric(GWS_signif) || length(GWS_signif) != n_period) {
    stop("GWS_signif must have length equal to length(period).", call. = FALSE)
  }
  if (!is.numeric(coi) || length(coi) != n_time) {
    stop("coi must have length equal to length(variable).", call. = FALSE)
  }

  # Ensure strictly increasing for boundary calculations
  if (any(!is.finite(period))) stop("'period' contains non-finite values.", call. = FALSE)
  if (is.unsorted(period, strictly = TRUE)) {
    stop("'period' must be strictly increasing.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Axes
  # ---------------------------------------------------------------------------
  time_year <- if (is.null(variable.year)) {
    seq_len(n_time)
  } else {
    if (!is.numeric(variable.year) || length(variable.year) != n_time) {
      stop("variable.year must be NULL or numeric with length equal to length(variable).", call. = FALSE)
    }
    if (is.unsorted(variable.year, strictly = TRUE)) {
      stop("variable.year must be strictly increasing.", call. = FALSE)
    }
    variable.year
  }

  # ---------------------------------------------------------------------------
  # Helper: compute cell boundaries for irregular grids
  # ---------------------------------------------------------------------------
  .cell_bounds <- function(x) {
    x <- as.numeric(x)
    mids <- (x[-1L] + x[-length(x)]) / 2
    lower <- c(x[1L] - (mids[1L] - x[1L]), mids)
    upper <- c(mids, x[length(x)] + (x[length(x)] - mids[length(mids)]))
    list(lower = lower, upper = upper)
  }

  xb <- .cell_bounds(time_year)
  yb <- .cell_bounds(period)

  # ---------------------------------------------------------------------------
  # Global Wavelet Spectrum data
  # ---------------------------------------------------------------------------
  GWS_df <- data.frame(
    period = as.numeric(period),
    GWS = as.numeric(GWS),
    GWS_signif = as.numeric(GWS_signif),
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Power field (native grid; explicit tile extents)
  # ---------------------------------------------------------------------------
  POWER_use <- as.matrix(POWER)
  POWER_use[!is.finite(POWER_use) | POWER_use <= 0] <- NA_real_
  power_log2 <- log2(POWER_use)

  # Build long dataframe with per-cell boundaries
  # Use t() so vectorization matches x-fastest ordering
  zvec <- as.vector(t(power_log2))  # length = n_time * n_period
  df_power <- data.frame(
    x = rep(time_year, times = n_period),
    y = rep(period, each = n_time),
    xmin = rep(xb$lower, times = n_period),
    xmax = rep(xb$upper, times = n_period),
    ymin = rep(yb$lower, each = n_time),
    ymax = rep(yb$upper, each = n_time),
    z = zvec,
    stringsAsFactors = FALSE
  )
  df_power <- df_power[is.finite(df_power$z), , drop = FALSE]
  if (nrow(df_power) == 0L) {
    stop("All POWER values are non-finite or <= 0 after cleaning; cannot plot.", call. = FALSE)
  }

  # Robust color limits
  zlims <- as.numeric(stats::quantile(df_power$z, probs = c(0.05, 0.95), na.rm = TRUE))
  if (!all(is.finite(zlims)) || zlims[1] >= zlims[2]) {
    zf <- range(df_power$z, finite = TRUE)
    zlims <- as.numeric(zf)
  }

  # ---------------------------------------------------------------------------
  # COI
  # ---------------------------------------------------------------------------
  coi_df <- data.frame(
    x = time_year,
    y = as.numeric(coi),
    stringsAsFactors = FALSE
  )
  pmin <- min(period, na.rm = TRUE)
  pmax <- max(period, na.rm = TRUE)
  coi_df <- coi_df[is.finite(coi_df$y) & (coi_df$y >= pmin) & (coi_df$y <= pmax), , drop = FALSE]

  # ---------------------------------------------------------------------------
  # Significance contour (native grid, no interpolation)
  # ---------------------------------------------------------------------------
  sigm_use <- as.matrix(sigm)
  sigm_use[!is.finite(sigm_use)] <- NA_real_

  sigm_df <- data.frame(
    x = rep(time_year, times = n_period),
    y = rep(period, each = n_time),
    z = as.vector(t(sigm_use)),
    stringsAsFactors = FALSE
  )
  sigm_df <- sigm_df[is.finite(sigm_df$z), , drop = FALSE]

  # ---------------------------------------------------------------------------
  # Panel a: Spectrum (continuous tiles)
  # ---------------------------------------------------------------------------
  p_spectrum <- ggplot2::ggplot(df_power) +
    ggplot2::theme_light() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = z),
      colour = NA
    ) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = zlims,
      oob = scales::squish
    ) +
    ggplot2::labs(x = "Time (year)", y = "Period (years)") +
    ggplot2::guides(fill = "none") +
    ggplot2::geom_line(
      data = coi_df,
      ggplot2::aes(x = x, y = y),
      linetype = "dashed",
      color = "red",
      linewidth = 0.85,
      inherit.aes = FALSE
    ) +
    ggplot2::stat_contour(
      data = sigm_df,
      ggplot2::aes(x = x, y = y, z = z),
      breaks = 1,
      color = "black",
      inherit.aes = FALSE
    )

  # ---------------------------------------------------------------------------
  # Panel b: Global Wavelet Spectrum
  # ---------------------------------------------------------------------------
  y_limits <- c(pmin, pmax)

  p_gws <- ggplot2::ggplot(GWS_df, ggplot2::aes(x = period, y = GWS)) +
    ggplot2::theme_light() +
    ggplot2::geom_line() +
    ggplot2::geom_point(shape = 19, size = 1) +
    ggplot2::geom_line(
      ggplot2::aes(x = period, y = GWS_signif),
      color = "red",
      linetype = "dashed",
      linewidth = 0.85
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_reverse(expand = c(0, 0), limits = y_limits) +
    ggplot2::coord_flip() +
    ggplot2::labs(y = bquote(Power ~ (.(variable.unit)^2)), x = "")

  # ---------------------------------------------------------------------------
  # Combine panels
  # ---------------------------------------------------------------------------
  p_spectrum + p_gws + patchwork::plot_layout(widths = c(3, 1.2))
}



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
