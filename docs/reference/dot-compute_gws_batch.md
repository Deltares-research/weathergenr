# Batch global wavelet spectrum for a simulated ensemble

Computes the unmasked global wavelet spectrum (GWS) for every column of
`sim_matrix` without allocating a three-dimensional wave array. Morlet
daughters are pre-computed once and reused across all series. Processing
is chunked to bound peak memory use.

The GWS formula matches `analyze_wavelet_spectrum(mode = "fast")`
exactly: each series is detrended (optional), standardized, zero-padded
to the next power of two, transformed via FFT, convolved with each
Morlet daughter via inverse FFT, and the unmasked GWS is
`variance_j * rowMeans(|wave_j|^2)`.

## Usage

``` r
.compute_gws_batch(
  sim_matrix,
  wavelet_pars,
  period = NULL,
  chunk_size = 5000L,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- sim_matrix:

  Numeric matrix. Rows are time steps, columns are realizations. Must
  not contain NA.

- wavelet_pars:

  Named list. Only `detrend` (logical) is used; the remaining entries
  (`signif_level`, `noise_type`, `period_lower_limit`) affect the
  observed series only and are ignored here.

- period:

  Numeric vector or NULL. Target period grid onto which each GWS is
  interpolated via
  [`stats::approx()`](https://rdrr.io/r/stats/approxfun.html). When
  NULL, or when the native grid already matches, no interpolation is
  performed. Provide the observed `period` vector from
  `analyze_wavelet_spectrum(obs_use)` here.

- chunk_size:

  Integer scalar. Columns per processing chunk. Smaller values reduce
  peak memory; larger values reduce loop overhead. Default 5000.

- parallel:

  Logical scalar. If TRUE, chunk-level GWS computation may run in
  parallel using a PSOCK cluster.

- n_cores:

  Integer scalar or NULL. Number of worker processes to use when
  `parallel = TRUE`. If NULL, all available logical cores may be used.

## Value

Numeric matrix of dimension `length(period) x ncol(sim_matrix)`
containing the unmasked GWS for each realization on `period`.
