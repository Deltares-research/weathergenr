#' Plot Wavelet Power Spectrum and Global Wavelet Spectrum
#'
#' Creates a two-panel figure showing the time-period wavelet power spectrum (with cone of influence and significance contour)
#' and the global wavelet spectrum for a time series.
#'
#' @param variable Numeric vector. The original time series (for length).
#' @param variable.year Numeric vector. Year labels for x-axis (optional).
#' @param period Numeric vector. Periods (scales) used in the wavelet transform.
#' @param POWER Matrix. Wavelet power spectrum (periods x time).
#' @param GWS Numeric vector. Global wavelet spectrum (mean power at each period).
#' @param GWS_signif Numeric vector. Significance threshold for the global wavelet spectrum.
#' @param coi Numeric vector. Cone of influence for each time point.
#' @param sigm Matrix. Significance mask (periods x time), e.g., significant (1) or not (0).
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

  # Input checks
  if (is.null(variable) || is.null(period) || is.null(POWER) ||
      is.null(GWS) || is.null(GWS_signif) || is.null(coi) || is.null(sigm)) {
    stop("All input arguments must be provided and non-null.")
  }

  n1 <- length(variable)
  n_period <- length(period)

  if (!is.matrix(POWER) || nrow(POWER) != n_period) {
    stop("POWER must be a matrix with nrow equal to length(period).")
  }
  if (!is.matrix(sigm) || nrow(sigm) != n_period) {
    stop("sigm must be a matrix with nrow equal to length(period).")
  }

  # Global Wavelet Spectrum data
  GWS_df <- data.frame(period = period, GWS = GWS, GWS_signif = GWS_signif,
     stringsAsFactors = FALSE)

  # Prepare time axis
  time_year <- if(is.null(variable.year)) {
    seq_len(length(variable))
  } else {
    variable.year
  }

  # Transform power to log scale
  power_log2 <- log(POWER, base = 2)
  power_log2_t <- t(power_log2)

  # Create data frame for power spectrum - METHOD 1: Manual construction
  n_time <- nrow(power_log2_t)
  n_per <- ncol(power_log2_t)

  df_power <- data.frame(
    x = rep(time_year, times = n_per),
    y = rep(period, each = n_time),
    z = as.vector(power_log2_t),
    stringsAsFactors = FALSE
  )

  # Interpolate for smoother contour
  interp <- suppressWarnings(akima::interp(
    x = df_power$x,
    y = df_power$y,
    z = df_power$z,
    xo = seq(min(df_power$x, na.rm = TRUE), max(df_power$x, na.rm = TRUE), length.out = 100),
    yo = seq(min(df_power$y, na.rm = TRUE), max(df_power$y, na.rm = TRUE), length.out = 100),
    extrap = TRUE,
    linear = FALSE
  ))

  # Convert interpolated data - CRITICAL FIX: Proper data frame construction
  interp_df <- data.frame(
    x = rep(interp$x, times = length(interp$y)),
    y = rep(interp$y, each = length(interp$x)),
    z = as.vector(interp$z),
    stringsAsFactors = FALSE
  )

  # Remove NAs that might cause issues
  interp_df <- interp_df[complete.cases(interp_df), ]

  # Cone of Influence (COI)
  coi_df <- data.frame(
    x = time_year,
    y = coi,
    stringsAsFactors = FALSE
  )
  coi_df <- coi_df[coi_df$y > min(interp_df$y, na.rm = TRUE), ]

  # Significance contour - METHOD 1: Manual construction
  sigm_t <- t(sigm)

  sigm_df <- data.frame(
    x = rep(time_year, times = n_per),
    y = rep(period, each = n_time),
    z = as.vector(sigm_t),
    stringsAsFactors = FALSE
  )

  # Common period scale
  pmin <- min(period, na.rm = TRUE)
  pmax <- max(period, na.rm = TRUE)
  y_limits <- c(pmin, pmax)

  # Panel a: Time-period power spectrum
  p_spectrum <- ggplot(interp_df, aes(x = x, y = y)) +
    theme_light() +
    geom_raster(aes(fill = z)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1) +
    labs(x = "Time (year)", y = "Period (years)") +
    guides(fill = "none") +
    geom_line(
      data = coi_df,
      aes(x = x, y = y),
      linetype = "dashed",
      color = "red",
      linewidth = 0.85,
      inherit.aes = FALSE
    ) +
    stat_contour(
      data = sigm_df,
      aes(x = x, y = y, z = z),
      breaks = c(-99, 1),
      color = "black",
      inherit.aes = FALSE
    )

  # Panel b: Global Wavelet Spectrum
  p_gws <- ggplot(GWS_df, aes(x = period, y = GWS)) +
    theme_light() +
    geom_line() +
    geom_point(shape = 19, size = 1) +
    geom_line(
      aes(x = period, y = GWS_signif),
      color = "red",
      linetype = "dashed",
      linewidth = 0.85
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_reverse(expand = c(0, 0), limits = y_limits) +
    coord_flip() +
    labs(y = bquote(Power ~ (.(variable.unit)^2)), x = "")

  # Combine panels
  p_combined <- p_spectrum + p_gws + patchwork::plot_layout(widths = c(3, 1.2))

  return(p_combined)
}




