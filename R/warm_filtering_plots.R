# ==============================================================================
# Internal Plotting Function for filter_warm_pool()
# ==============================================================================

#' Create diagnostic plots for filtered pool
#'
#' @description
#' Internal function used by filter_warm_pool() to create diagnostic plots.
#' Shows time series, statistics scatter, and wavelet GWS comparisons.
#'
#' @param obs_series Numeric vector of observed values
#' @param sim_series Numeric matrix of simulated values
#' @param pool Integer vector of pool indices
#' @param rel_diff_mean Relative differences in mean
#' @param rel_diff_sd Relative differences in SD
#' @param tail_metrics Tail metrics list
#' @param power_period Wavelet periods
#' @param power_obs Observed GWS
#' @param power_signif Significance curve
#' @param gws_cache Cached GWS matrix (n_periods x n_realizations)
#' @param wavelet_q Two quantiles for ribbon (e.g., c(0.50, 0.95))
#'
#' @return List of ggplot objects (timeseries, stats, wavelet_gws)
#'
#' @keywords internal
#' @import ggplot2
#' @importFrom stats quantile
#' @export
plot_filter_diagnostics <- function(obs_series, sim_series, pool,
                                    rel_diff_mean, rel_diff_sd, tail_metrics,
                                    power_period, power_obs, power_signif,
                                    gws_cache, wavelet_q = c(0.05, 0.95)) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' required for plotting.", call. = FALSE)
  }

  # Validate inputs
  if (!is.numeric(obs_series) || !is.vector(obs_series)) {
    stop("obs_series must be numeric vector.", call. = FALSE)
  }
  if (!is.matrix(sim_series) || !is.numeric(sim_series)) {
    stop("sim_series must be numeric matrix.", call. = FALSE)
  }
  if (nrow(sim_series) != length(obs_series)) {
    stop("sim_series rows must match length(obs_series).", call. = FALSE)
  }
  if (!is.numeric(wavelet_q) || length(wavelet_q) != 2L ||
      any(!is.finite(wavelet_q)) || any(wavelet_q <= 0) || any(wavelet_q >= 1) ||
      wavelet_q[1] >= wavelet_q[2]) {
    stop("wavelet_q must be two probabilities in (0,1) with q1<q2.", call. = FALSE)
  }

  n_pool <- length(pool)

  # Extract pool data
  sim_pool <- sim_series[, pool, drop = FALSE]
  if (is.null(colnames(sim_pool))) colnames(sim_pool) <- paste0("rlz_", pool)

  # Compute approximate relative differences for plotting
  # Convert log-ratio back to percentage for readability
  e <- tail_metrics$M_obs_low  # Use small epsilon
  if (!is.finite(e)) e <- 1e-5

  lr_low  <- log((tail_metrics$M_sim_low[pool]  + e) / (tail_metrics$M_obs_low  + e))
  lr_high <- log((tail_metrics$M_sim_high[pool] + e) / (tail_metrics$M_obs_high + e))

  stats_df <- data.frame(idx = pool, mean = rel_diff_mean[pool] * 100,
                         sd   = rel_diff_sd[pool]   * 100, tail_low  = (exp(lr_low)  - 1) * 100,
                         tail_high = (exp(lr_high) - 1) * 100)

  stats_pool_long <- stats_df |>
    tidyr::pivot_longer(cols = -idx, names_to = "par", values_to = "value")

  # Extract pool GWS from cache
  gws_pool_mat <- gws_cache[, pool, drop = FALSE]
  gws_pool_mean <- fill_nearest(rowMeans(gws_pool_mat, na.rm = TRUE))
  gws_pool_q_lo <- fill_nearest(apply(gws_pool_mat, 1, stats::quantile,
                                      probs = wavelet_q[1], na.rm = TRUE, names = FALSE, type = 4))
  gws_pool_q_hi <- fill_nearest(apply(gws_pool_mat, 1, stats::quantile,
                                      probs = wavelet_q[2], na.rm = TRUE, names = FALSE, type = 4))

  # ===========================================================================
  # Time series plot
  # ===========================================================================

  df_ts_sim_long <- data.frame(
    year = seq_len(length(obs_series)),
    as.data.frame(sim_pool, check.names = FALSE),
    check.names = FALSE) |>
    tidyr::pivot_longer(cols = -year, names_to = "series", values_to = "value")

  df_ts_obs <- data.frame(year = seq_len(length(obs_series)), value = obs_series)

  p_timeseries <- ggplot() + theme_light() +
    geom_line(data = df_ts_sim_long, aes(x = year, y = value, group = series),
              linewidth = 0.7, alpha = 0.20) +
    geom_line(data = df_ts_obs, aes(x = year, y = value),
              color = "blue", linewidth = 0.9) +
    labs(x = "Year index", y = "Value")

  # ===========================================================================
  # Statistics scatter plot
  # ===========================================================================
  par_labels <- c(mean = "Mean", sd = "Standard\nDeviation",
                  tail_low = "Dry extremes\n(lower-tail mass)", tail_high = "Wet extremes\n(upper-tail mass)")

  stats_pool_long <- stats_df |>
    tidyr::pivot_longer(cols = -idx, names_to = "par", values_to = "value") |>
    dplyr::filter(is.finite(.data$value))

  stats_plot <- stats_pool_long |>
    dplyr::mutate(par = factor(.data$par, levels = c("mean", "sd", "tail_low", "tail_high")))


  p_stats <- ggplot(stats_plot, aes(x = "a", y = .data$value)) +
    theme_bw() +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "blue") +
    geom_point(size = 3, alpha = 0.3, color = "black", position = "jitter") +
    facet_wrap(~par, nrow = 1, labeller = as_labeller(par_labels)) +
    labs(x = "Filtered realizations", y = "Relative difference (%)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(-50, 50), breaks = seq(-50, 50, 25))

  # ===========================================================================
  # Wavelet GWS plot
  # ===========================================================================

  df_gws_sum <- data.frame(period = power_period, mean = gws_pool_mean, lo = gws_pool_q_lo,
                           hi = gws_pool_q_hi)

  df_gws_lines <- data.frame(period = power_period, obs = as.numeric(power_obs),
                             signif = as.numeric(power_signif), curve = "Significance")

  p_gws <- ggplot() + theme_light() +
    geom_ribbon(aes(x = period, ymin = lo, ymax = hi), df_gws_sum, alpha = 0.20) +
    geom_line(aes(x = period, y = mean), df_gws_sum, linewidth = 0.8) +
    geom_line(aes(x = period, y = obs), df_gws_lines, color = "blue", linewidth = 0.8) +
    geom_line(aes(x = period, y = signif), df_gws_lines, color = "red", linewidth = 0.8, linetype = "dashed") +
    labs(x = "Period", y = "Power")

  list(timeseries = p_timeseries, stats = p_stats, wavelet_gws = p_gws)
}
