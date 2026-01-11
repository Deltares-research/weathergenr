# Compute tail-mass metrics for filtering

Computes tail-mass metrics based on lower and upper tail quantile
thresholds. Uses robust scale estimation (IQR -\> MAD -\> SD fallback)
and normalizes tail deficit/excess masses by series length and scale.

## Usage

``` r
compute_tailmass_metrics(
  obs.use,
  series_sim_for_stats,
  tail.low.p,
  tail.high.p,
  tail.eps
)
```

## Arguments

- obs.use:

  Numeric vector of observed values

- series_sim_for_stats:

  Numeric matrix of simulated values (n_use x n_realizations)

- tail.low.p:

  Lower tail quantile probability (e.g., 0.10)

- tail.high.p:

  Upper tail quantile probability (e.g., 0.90)

- tail.eps:

  Epsilon for log transform to avoid log(0)

## Value

List with elements:

- thr_low, thr_high: Observed quantile thresholds

- scale_obs: Robust scale estimate

- M_obs_low, M_obs_high: Observed tail masses

- M_sim_low, M_sim_high: Simulated tail masses (vectors)

- logdiff_low, logdiff_high: Log-distance metrics (vectors)
