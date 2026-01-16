# Compute tail-mass metrics for filtering

Computes tail-mass metrics based on lower and upper tail quantile
thresholds. Uses robust scale estimation (IQR -\> MAD -\> SD fallback)
and normalizes tail deficit/excess masses by series length and scale.

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

  Numeric vector of observed values.

- sim_series_stats:

  Numeric matrix of simulated values (n_use x n_realizations).

- tail_low_p:

  Lower tail quantile probability (e.g., 0.10).

- tail_high_p:

  Upper tail quantile probability (e.g., 0.90).

- tail_eps:

  Epsilon for log transform to avoid log(0).

## Value

List with tail thresholds, scale, masses, and log-distance metrics.
