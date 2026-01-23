# Filter and sample WARM realizations using distributional, tail, and spectral criteria

Filters an ensemble of annual WARM realizations against an observed
annual reference series. Filtering is based on three families of
criteria: distributional moments, tail mass behavior, and spectral
similarity from the global wavelet spectrum.

The function computes pass or fail vectors for each filter family and
builds a candidate pool. If the pool is smaller than the requested
sample size, the function relaxes thresholds iteratively until enough
candidates are found or the maximum number of relaxation iterations is
reached.

Relaxation is adaptive. At each iteration, the currently most
restrictive active filter is relaxed, defined as the filter with the
lowest pass rate among the active filters. This is constrained by the
relaxation order given in `relax_order`, which sets which filters are
eligible to be relaxed and how wavelet relaxation is parameterized.

If the pool is still smaller than `n_select` after relaxation, a
deterministic fallback selects the `n_select` realizations with the
smallest absolute relative mean difference.

## Usage

``` r
filter_warm_pool(
  obs_series = NULL,
  sim_series = NULL,
  n_select = 5,
  seed = NULL,
  relax_order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
  filter_bounds = list(),
  wavelet_args = list(signif_level = 0.8, noise_type = "red", period_lower_limit = 2,
    detrend = TRUE),
  modwt_n_levels = NULL,
  make_plots = FALSE,
  parallel = FALSE,
  n_cores = NULL,
  cache_gws = FALSE,
  verbose = FALSE
)
```

## Arguments

- obs_series:

  Numeric vector. Observed annual series used as the reference target.

- sim_series:

  Numeric matrix. Simulated annual realizations with years in rows and
  realizations in columns.

- n_select:

  Integer scalar. Number of realizations to return in the `selected`
  element.

- seed:

  Integer scalar or NULL. If provided, sets the random seed used for
  selecting alignment windows when `obs_series` and `sim_series` have
  different lengths, and for sampling within the final candidate pool.

- relax_order:

  Character vector. Relaxation priority ordering for the filter
  families. Must contain exactly: `"mean"`, `"sd"`, `"tail_low"`,
  `"tail_high"`, `"wavelet"`.

- filter_bounds:

  Named list. Overrides for filtering thresholds and relaxation
  controls. Keys must be snake case and match those returned by
  [`filter_warm_bounds_defaults()`](https://deltares-research.github.io/weathergenr/reference/filter_warm_bounds_defaults.md).

- wavelet_args:

  Named list. Parameters passed to
  [`analyze_wavelet_spectrum()`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md)
  for observed and simulated series. Expected entries include
  `signif_level`, `noise_type`, `period_lower_limit`, and `detrend`.

- modwt_n_levels:

  Integer or NULL. Number of MODWT levels used in WARM. This value is
  used only for diagnostics. If NULL, a value is estimated from the
  series length.

- make_plots:

  Logical scalar. If TRUE, compute diagnostic plots for the selected
  realizations and return them in `plots`.

- parallel:

  Logical scalar. If TRUE, allows spectral metrics to be computed in
  parallel inside
  [`compute_spectral_metrics()`](https://deltares-research.github.io/weathergenr/reference/compute_spectral_metrics.md).

- n_cores:

  Integer scalar or NULL. Number of worker processes to use for parallel
  execution. If NULL, an internal default is used.

- cache_gws:

  Logical scalar. If TRUE, cache simulated global wavelet spectra for
  diagnostics; this can increase memory use.

- verbose:

  Logical scalar. If TRUE, logs setup information, per iteration pass
  rates, and relaxation actions.

## Value

Named list with the following elements:

- pool:

  Numeric matrix. Candidate pool of realizations that passed the final
  set of filters. Subset of columns from `sim_series`.

- selected:

  Numeric matrix. The `n_select` realizations chosen from the pool.
  Subset of columns from `sim_series`.

- summary:

  Data frame. Pass counts and pass rates for each filter family, plus
  the selection mode used.

- diagnostics:

  List. Window metadata, relaxation log, spectral diagnostics, spectral
  metrics, and the final bounds used. When `cache_gws = TRUE`, includes
  `gws_cache`.

- plots:

  NULL or a named list of ggplot objects when `make_plots = TRUE`.
