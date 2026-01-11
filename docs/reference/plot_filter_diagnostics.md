# Create diagnostic plots for filtered pool

Internal function used by filter_warm_simulations() to create diagnostic
plots. Shows time series, statistics scatter, and wavelet GWS
comparisons.

## Usage

``` r
plot_filter_diagnostics(
  obs.use,
  series_sim_for_stats,
  pool,
  rel.diff.mean,
  rel.diff.sd,
  tail_metrics,
  power.period,
  power.obs,
  power.signif,
  gws_cache_mat,
  wavelet_q = c(0.05, 0.95)
)
```

## Arguments

- obs.use:

  Numeric vector of observed values

- series_sim_for_stats:

  Numeric matrix of simulated values

- pool:

  Integer vector of pool indices

- rel.diff.mean:

  Relative differences in mean

- rel.diff.sd:

  Relative differences in SD

- tail_metrics:

  Tail metrics list

- power.period:

  Wavelet periods

- power.obs:

  Observed GWS

- power.signif:

  Significance curve

- gws_cache_mat:

  Cached GWS matrix (n_periods x n_realizations)

- wavelet_q:

  Two quantiles for ribbon (e.g., c(0.50, 0.95))

## Value

List of ggplot objects (timeseries, stats, wavelet_gws)
