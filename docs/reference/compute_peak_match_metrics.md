# Compute significant peak match metrics for one simulated spectrum

Compares a simulated global wavelet spectrum to a set of significant
observed peaks. For each observed peak, a match is declared if the
simulated spectrum contains a local maximum within a period tolerance
and the matched power is within a magnitude tolerance.

Period tolerance is measured in log2 period space. Magnitude tolerance
is measured as the absolute log ratio between simulated and observed
power.

## Usage

``` r
compute_peak_match_metrics(
  gws_sim,
  period,
  obs_peaks,
  period_tol = 0.5,
  mag_tol_log = log(1.5),
  eps = 1e-10,
  log2_period = NULL
)
```

## Arguments

- gws_sim:

  Numeric vector. Simulated global wavelet spectrum aligned to `period`.

- period:

  Numeric vector. Period grid for `gws_sim`.

- obs_peaks:

  Data frame. Significant observed peaks returned by
  [`identify_significant_peaks()`](https://deltares-research.github.io/weathergenr/reference/identify_significant_peaks.md).

- period_tol:

  Numeric scalar. Maximum allowed absolute difference in log2 period
  between an observed peak period and a simulated candidate period.

- mag_tol_log:

  Numeric scalar. Maximum allowed absolute log ratio between simulated
  and observed peak power.

- eps:

  Numeric scalar. Positive constant used for numerical stability.

- log2_period:

  Numeric vector or NULL. Pre-computed log2(period) for efficiency. If
  NULL, computed internally.

## Value

Named list with: `peak_match_frac` fraction of observed peaks that were
matched, `peak_mag_mean_abs_log_ratio` mean absolute log ratio for
matched peaks, or NA if no peaks were matched.
