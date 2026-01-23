# Compute tail mass metrics for filtering

Computes tail mass metrics for an observed series and a simulated
ensemble. Tail mass metrics quantify how much probability mass lies
beyond observed quantile thresholds, expressed as normalized deficit or
excess mass.

The method uses robust scale estimation for normalization. The scale is
chosen as IQR when available, otherwise MAD, otherwise standard
deviation, otherwise 1. Tail masses are normalized by series length and
scale to make values more comparable across datasets.

The returned log difference vectors are used by
[`filter_warm_pool()`](https://deltares-research.github.io/weathergenr/reference/filter_warm_pool.md)
to decide whether a realization passes tail filters.

## Usage

``` r
compute_tailmass_metrics(
  obs_use,
  sim_series_stats,
  tail_low_p,
  tail_high_p,
  tail_eps
)
```

## Arguments

- obs_use:

  Numeric vector. Observed annual values after any window alignment
  performed by the caller.

- sim_series_stats:

  Numeric matrix. Simulated values aligned to `obs_use`. Rows are years
  and columns are realizations.

- tail_low_p:

  Numeric scalar. Lower tail quantile probability used to compute the
  low threshold on the observed series.

- tail_high_p:

  Numeric scalar. Upper tail quantile probability used to compute the
  high threshold on the observed series.

- tail_eps:

  Numeric scalar. Positive constant added inside log transforms to avoid
  log of zero.

## Value

Named list with tail thresholds, scale, observed and simulated tail mass
metrics, and log difference vectors:

- thr_low:

  Lower threshold from the observed series.

- thr_high:

  Upper threshold from the observed series.

- scale_obs:

  Robust scale used for normalization.

- M_obs_low:

  Observed low tail deficit mass.

- M_obs_high:

  Observed high tail excess mass.

- M_sim_low:

  Vector of simulated low tail deficit masses.

- M_sim_high:

  Vector of simulated high tail excess masses.

- logdiff_low:

  Vector of absolute log differences for low tail mass.

- logdiff_high:

  Vector of absolute log differences for high tail mass.
