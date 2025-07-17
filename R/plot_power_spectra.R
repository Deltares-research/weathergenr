#' Plot Wavelet Power Spectrum and Global Wavelet Spectrum
#'
#' Creates a two-panel figure showing the time-period wavelet power spectrum (with cone of influence and significance contour)
#' and the global wavelet spectrum for a time series.
#'
#' @param variable Numeric vector. The original time series (for length).
#' @param period Numeric vector. Periods (scales) used in the wavelet transform.
#' @param POWER Matrix. Wavelet power spectrum (periods x time).
#' @param GWS Numeric vector. Global wavelet spectrum (mean power at each period).
#' @param GWS_signif Numeric vector. Significance threshold for the global wavelet spectrum.
#' @param coi Numeric vector. Cone of influence for each time point.
#' @param sigm Matrix. Significance mask (periods x time), e.g., significant (1) or not (0).
#' @param variable.unit Character. Unit label for the variable (default is "mm").
#'
#' @import ggplot2
#' @return A patchwork plot object with the time-period power spectrum and global wavelet spectrum.
#' @export
#'
#' @examples
#' # plot_wavelet_spectra(variable, period, POWER, GWS, GWS_signif, coi, sigm, variable.unit = "mm")
plot_wavelet_spectra <- function(
    variable,
    period,
    POWER,
    GWS,
    GWS_signif,
    coi,
    sigm,
    variable.unit = "mm"
) {

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
  GWS_df <- data.frame(period = period, GWS = GWS, GWS_signif = GWS_signif)

  # Prepare data for power spectrum plot (Panel a)
  power_log2 <- t(log(POWER, base = 2))
  df_power <- tibble::as_tibble(power_log2, .name_repair = ~ as.character(period)) %>%
    dplyr::mutate(x = 1:length(variable)) %>%
    tidyr::pivot_longer(-x, names_to = "y", values_to = "z") %>%
    dplyr::mutate(y = as.numeric(y), z = as.numeric(z), x = as.numeric(x))

  # Interpolate for smoother contour (optional, can skip if not needed)
  interp <- suppressWarnings(akima::interp(
    x = df_power$x, y = df_power$y, z = df_power$z,
    xo = seq(min(df_power$x), max(df_power$x), length = 20),
    yo = seq(min(df_power$y), max(df_power$y), length = 20),
    extrap = TRUE, linear = FALSE))

  interp_df <- tibble::as_tibble(interp$z, .name_repair = ~ as.character(interp$y)) %>%
    dplyr::mutate(x = interp$x) %>%
    tidyr::pivot_longer(-x, names_to = "y", values_to = "z") %>%
    dplyr::mutate(y = as.numeric(y), z = as.numeric(z), x = as.numeric(x))

  # Cone of Influence (COI)
  coi_df <- tibble::tibble(x = 1:n1, y = coi) %>%
    dplyr::filter(y > min(interp_df$y, na.rm = TRUE))

  # Significance contour (Panel a, black contour)
  sigm_t <- t(sigm)
  sigm_df <- tibble::as_tibble(sigm_t, .name_repair = ~ as.character(period)) %>%
    dplyr::mutate(x = 1:n1) %>%
    tidyr::pivot_longer(-x, names_to = "y", values_to = "z") %>%
    dplyr::mutate(y = as.numeric(y), z = as.numeric(z), x = as.numeric(x))

  # Panel a: Time-period power spectrum
  p_spectrum <- ggplot(interp_df, aes(x = x, y = y)) +
    theme_light() +
    geom_raster(aes(fill = z)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fill_distiller(palette = "YlGnBu") +
    labs(x = "Time (years)", y = "Period (years)") +
    guides(fill = "none") +
    geom_line(
      data = coi_df, aes(x = x, y = y),
      linetype = "dashed", color = "red", linewidth = 0.85
    ) +
    stat_contour(
      data = sigm_df,
      aes(z = z, x = x, y = y),
      breaks = c(-99, 1), color = "black"
    )

  # Panel b: Global Wavelet Spectrum
  p_gws <- ggplot(GWS_df) +
    theme_light() +
    geom_line(aes(period, GWS)) +
    geom_point(aes(period, GWS), shape = 19, size = 1) +
    geom_line(aes(period, GWS_signif),
                       color = "red", linetype = "dashed", linewidth = 0.85) +
    scale_y_continuous() +
    scale_x_reverse(expand = c(0, 0)) +
    coord_flip() +
    labs(y = bquote(Power ~ (.(variable.unit)^2)), x = "")

  # Combine panels
  p_combined <- p_spectrum + p_gws + patchwork::plot_layout(widths = c(2, 1.25))

  return(p_combined)
}

