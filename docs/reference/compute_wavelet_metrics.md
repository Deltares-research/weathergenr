# Compute wavelet metrics for all realizations

Performs wavelet analysis on observed series and all simulated
realizations. Identifies significant periods, computes regional power
and contrast metrics, and ALWAYS caches both masked (for filtering) and
unmasked (for plotting) GWS.

## Usage

``` r
compute_wavelet_metrics(
  obs.use,
  series_sim_for_stats,
  wavelet.pars,
  padding,
  min_bg
)
```

## Arguments

- obs.use:

  Numeric vector of observed values

- series_sim_for_stats:

  Numeric matrix of simulated values

- wavelet.pars:

  List of wavelet parameters (signif.level, noise.type, etc.)

- padding:

  Logical for period padding

- min_bg:

  Minimum background power threshold

## Value

List with elements:

- active: Logical, whether wavelet filter is active

- diagnostics: Detailed wavelet diagnostics

- power.period, power.obs, power.signif: Period and power vectors

- P_sim_reg: Regional power matrix (n_realizations x n_regions)

- P_sim_bg: Background power vector (n_realizations)

- presence_rpad: Presence indicators (n_realizations)

- gws_cache: ALWAYS cached masked GWS matrix (n_periods x
  n_realizations) for filtering

- gws_cache_unmasked: ALWAYS cached unmasked GWS matrix (n_periods x
  n_realizations) for plotting
