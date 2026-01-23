# Compute spectral match metrics for one realization

Computes two spectral similarity metrics between an observed and a
simulated global wavelet spectrum: 1) correlation between log
transformed spectra, and 2) significant peak match metrics relative to
significant observed peaks.

## Usage

``` r
compute_spectral_match_single(
  gws_obs,
  gws_sim,
  period,
  obs_peaks,
  peak_period_tol = 0.5,
  peak_mag_tol_log = log(1.5),
  eps = 1e-10,
  log2_period = NULL
)
```

## Arguments

- gws_obs:

  Numeric vector. Observed global wavelet spectrum aligned to `period`.

- gws_sim:

  Numeric vector. Simulated global wavelet spectrum aligned to `period`.

- period:

  Numeric vector. Period grid for both spectra.

- obs_peaks:

  Data frame. Significant observed peaks returned by
  [`identify_significant_peaks()`](https://deltares-research.github.io/weathergenr/reference/identify_significant_peaks.md).

- peak_period_tol:

  Numeric scalar. Period matching tolerance in log2 period space.

- peak_mag_tol_log:

  Numeric scalar. Magnitude matching tolerance as an absolute log ratio.

- eps:

  Numeric scalar. Positive constant used to avoid log of zero.

- log2_period:

  Numeric vector or NULL. Pre-computed log2(period) for efficiency. If
  NULL, computed internally.

## Value

Named list with: `spectral_cor` correlation between log transformed
spectra, `peak_match_frac` fraction of significant observed peaks
matched, `peak_mag_mean_abs_log_ratio` mean absolute log ratio for
matched peaks.
