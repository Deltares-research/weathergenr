# Compute spectral metrics for an ensemble of realizations

Computes wavelet based spectral similarity metrics between an observed
annual series and each realization in a simulated ensemble.

The observed series is analyzed once to obtain a period grid, an
observed global wavelet spectrum, a significance curve, and a set of
significant peaks. Each simulated realization is then analyzed on the
same period grid and is scored using: 1) correlation between log
transformed spectra, and 2) significant peak match fraction and
magnitude agreement.

If `parallel = TRUE`, the per realization computation can run in
parallel. Parallel execution uses a PSOCK cluster created by the base
parallel package, which is typically the most portable option across
operating systems.

## Usage

``` r
compute_spectral_metrics(
  obs_use,
  sim_series_stats,
  wavelet_pars,
  modwt_n_levels = NULL,
  n_sig_peaks_max = 2L,
  peak_period_tol = 0.5,
  peak_mag_tol_log = log(1.5),
  eps = 1e-10,
  parallel = FALSE,
  n_cores = NULL,
  cache_gws = FALSE
)
```

## Arguments

- obs_use:

  Numeric vector. Observed annual series after any window alignment
  performed by the caller.

- sim_series_stats:

  Numeric matrix. Simulated ensemble after any window alignment
  performed by the caller. Rows are years and columns are realizations.

- wavelet_pars:

  Named list. Parameters passed to
  [`analyze_wavelet_spectrum()`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md).
  Expected entries include `signif_level`, `noise_type`,
  `period_lower_limit`, and `detrend`.

- modwt_n_levels:

  Integer or NULL. Number of MODWT levels used in WARM. Used only for
  diagnostics. If NULL, a value is estimated from series length.

- n_sig_peaks_max:

  Integer scalar. Maximum number of significant observed peaks to
  enforce when computing peak matching metrics.

- peak_period_tol:

  Numeric scalar. Period matching tolerance in log2 period space.

- peak_mag_tol_log:

  Numeric scalar. Magnitude matching tolerance as an absolute log ratio.

- eps:

  Numeric scalar. Positive constant used for numerical stability in log
  transforms and ratios.

- parallel:

  Logical scalar. If TRUE, compute per realization spectra and metrics
  in parallel.

- n_cores:

  Integer scalar or NULL. Number of worker processes to use when
  `parallel = TRUE`. If NULL, an internal default is used.

- cache_gws:

  Logical scalar. If TRUE, store simulated global wavelet spectra in
  `gws_cache`. When FALSE, `gws_cache` is returned as NULL and simulated
  spectra are not retained.

## Value

Named list with:

- active:

  Logical scalar. TRUE when wavelet metrics were computed.

- period:

  Numeric vector. Period grid used for all spectra.

- gws_obs:

  Numeric vector. Observed global wavelet spectrum on `period`.

- gws_signif:

  Numeric vector or NULL. Significance curve on `period`.

- gws_cache:

  Numeric matrix or NULL. Simulated spectra on `period` for each
  realization when `cache_gws = TRUE`; otherwise NULL.

- metrics:

  List with numeric vectors `spectral_cor`, `peak_match_frac`, and
  `peak_mag_mean_abs_log_ratio`.

- diagnostics:

  List. Summary information including significant peaks and observed
  spectrum summary statistics.
