# Compute wavelet metrics for all realizations

performs wavelet analysis on observed series and all simulated
realizations. Identifies significant periods, computes regional power
and contrast metrics, and caches both masked (for filtering) and
unmasked (for plotting) GWS.

## Usage

``` r
compute_wavelet_metrics(
  obs_use,
  sim_series_stats,
  wavelet_pars,
  padding,
  min_bg
)
```

## Arguments

- obs_use:

  Numeric vector of observed values.

- sim_series_stats:

  Numeric matrix of simulated values.

- wavelet_pars:

  List of wavelet parameters (signif_level, noise_type, etc.).

- padding:

  Logical for period padding.

- min_bg:

  Minimum background power threshold.

## Value

List with wavelet filter diagnostics and cached spectra.
