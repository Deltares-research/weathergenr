# Create diagnostic plots for filtered pool

Internal function used by filter_warm_pool() to create diagnostic plots.
Shows time series, statistics scatter, and wavelet GWS comparisons.

## Usage

``` r
plot_filter_diagnostics(
  obs_series,
  sim_series,
  pool,
  rel_diff_mean,
  rel_diff_sd,
  tail_metrics,
  power_period,
  power_obs,
  power_signif,
  wavelet_pars,
  wavelet_q = c(0.05, 0.95)
)
```

## Arguments

- obs_series:

  Numeric vector of observed values

- sim_series:

  Numeric matrix of simulated values

- pool:

  Integer vector of pool indices

- rel_diff_mean:

  Relative differences in mean

- rel_diff_sd:

  Relative differences in SD

- tail_metrics:

  Tail metrics list

- power_period:

  Wavelet periods

- power_obs:

  Observed GWS

- power_signif:

  Significance curve

- wavelet_pars:

  Named list. Wavelet settings passed to analyze_wavelet_spectrum()
  (signif_level, noise_type, period_lower_limit, detrend).

- wavelet_q:

  Two quantiles for ribbon (e.g., c(0.50, 0.95))

## Value

List of ggplot objects (timeseries, stats, wavelet_gws)
