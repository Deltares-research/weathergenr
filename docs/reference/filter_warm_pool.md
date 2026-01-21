# Filter and sample WARM realizations using distributional, tail, and wavelet criteria

Filters an ensemble of WARM-generated annual realizations against an
observed annual series using three criterion families:

- **Distributional**: relative differences in mean and standard
  deviation.

- **Tail behaviour**: lower/upper tail mass relative to observed
  quantile thresholds.

- **Spectral**: observed-relevant global wavelet spectrum (GWS)
  filtering.

If fewer than `n_select` realizations pass, the function relaxes
criteria iteratively (up to `filter_bounds$relax_max_iter`) by loosening
the currently most restrictive active filter (lowest pass rate). If
still insufficient, a deterministic fallback returns the `n_select`
realizations with the smallest absolute relative mean difference.

## Usage

``` r
filter_warm_pool(
  obs_series = NULL,
  sim_series = NULL,
  n_select = 5,
  seed = NULL,
  pad_periods = TRUE,
  relax_order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
  filter_bounds = list(),
  wavelet_args = list(signif_level = 0.8, noise_type = "red", period_lower_limit = 2,
    detrend = TRUE),
  make_plots = FALSE,
  verbose = FALSE
)
```

## Arguments

- obs_series:

  Numeric vector. Observed annual series used as the reference.

- sim_series:

  Numeric matrix. Simulated annual realizations with years in rows and
  realizations in columns.

- n_select:

  Integer scalar. Number of realizations to return in `selected`.

- seed:

  Optional integer scalar. Random seed used for window selection (if
  lengths differ) and for sampling from the final candidate pool.

- pad_periods:

  Logical scalar. If `TRUE`, expands the observed significant-period
  band by one index on each side when checking simulated presence in the
  observed-relevant band.

- relax_order:

  Character vector. Relaxation priority ordering for criteria. Must
  contain each of `c("mean","sd","tail_low","tail_high","wavelet")`
  exactly once.

- filter_bounds:

  Named list. Filtering thresholds and relaxation controls. Any entry
  overrides internal defaults. Uses snake_case keys (e.g. `tail_low_p`,
  not `tail.low.p`).

- wavelet_args:

  Named list passed to
  [`analyze_wavelet_spectrum`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md)
  for observed and simulated series (e.g., `signif_level`, `noise_type`,
  `period_lower_limit`, `detrend`).

- make_plots:

  Logical scalar. If `TRUE`, returns diagnostic plots in `plots`.

- verbose:

  Logical scalar. If `TRUE`, logs per-iteration pass rates and
  relaxation steps.

## Value

A list with:

- pool:

  Numeric matrix. Final candidate pool (subset of columns from
  `sim_series`).

- selected:

  Numeric matrix. `n_select` realizations selected from the final pool.

- summary:

  Data frame summarising pass counts/rates and selection mode.

- diagnostics:

  List with window metadata, indices, relaxation log, and final bounds.

- plots:

  NULL or a named list of ggplot objects when `make_plots = TRUE`.
