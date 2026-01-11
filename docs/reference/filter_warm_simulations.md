# Filter and sample WARM simulations using distributional, tail, and wavelet criteria

Filters a set of WARM realizations (simulated annual series) against an
observed annual series using three tiers of checks:

- **Distributional**: relative differences in mean and standard
  deviation,

- **Tail behavior**: lower/upper tail mass relative to observed quantile
  thresholds,

- **Spectral**: an observed-relevant global wavelet spectrum (GWS)
  filter.

If fewer than `sample.num` realizations pass, the function relaxes
criteria iteratively (up to `bounds$relax.max.iter`) by loosening the
currently most restrictive active filter. If still insufficient, a
deterministic fallback returns the `sample.num` realizations with the
smallest absolute relative mean difference.

## Usage

``` r
filter_warm_simulations(
  series.obs = NULL,
  series.sim = NULL,
  sample.num = 5,
  seed = NULL,
  padding = TRUE,
  relax.order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
  bounds = list(),
  wavelet.pars = list(signif.level = 0.8, noise.type = "red", period.lower.limit = 2,
    detrend = TRUE),
  make.plots = FALSE,
  verbose = FALSE
)
```

## Arguments

- series.obs:

  Numeric vector. Observed annual series (no missing values
  recommended).

- series.sim:

  Numeric matrix. Simulated annual series with years in rows and
  realizations in columns.

- sample.num:

  Integer scalar. Number of realizations to return in `sampled`. If
  greater than `ncol(series.sim)`, it is reduced to `ncol(series.sim)`
  with a warning.

- seed:

  Optional integer. Random seed used for window selection (if lengths
  differ) and for sampling from the final candidate pool.

- padding:

  Logical scalar. If `TRUE`, expands the observed significant-period set
  by one index on each side when assessing presence of observed-relevant
  signal in simulated spectra.

- relax.order:

  Character vector. Relaxation priority ordering for criteria. Must
  contain each of `c("mean","sd","tail_low","tail_high","wavelet")`
  exactly once. When multiple filters have identical pass rates, this
  order is used to break ties.

- bounds:

  Named list. Filtering thresholds and relaxation controls. Any entry
  provided overrides the defaults. Common entries include:

  - `mean`, `sd`: relative-difference tolerances for mean and sd.

  - `tail.low.p`, `tail.high.p`: quantiles defining lower/upper tail
    thresholds.

  - `tail.tol.log`: tolerance on log-distance between simulated and
    observed tail mass.

  - `tail.eps`: epsilon added to tail mass before log transform for
    stability.

  - `sig.frac`: minimum fraction of observed-relevant regions required
    to pass the wavelet filter.

  - `wavelet.region.tol`, `wavelet.contrast.tol`: tolerances for
    regional power and contrast ratios.

  - `wavelet.require_presence`, `wavelet.presence.frac`: controls for
    requiring signal presence.

  - `relax.*`: relaxation step sizes and caps (including
    `relax.max.iter`).

- wavelet.pars:

  Named list passed to
  [`wavelet_spectral_analysis`](https://github.com/Deltares-research/weathergenr/reference/wavelet_spectral_analysis.md)
  for observed and simulated series. Typical entries include
  `signif.level`, `noise.type`, `period.lower.limit`, and `detrend`.

- make.plots:

  Logical scalar. If `TRUE`, returns diagnostic plots (time series
  overlay, relative-difference summaries, and wavelet GWS diagnostics)
  in `plots`.

- verbose:

  Logical scalar. If `TRUE`, logs per-iteration pass rates and
  relaxation steps using
  [`logger::log_info()`](https://daroczig.github.io/logger/reference/log_level.html).

## Value

A list with:

- subsetted:

  Numeric matrix. All realizations that remain in the final candidate
  pool (columns subset of `series.sim`).

- sampled:

  Numeric matrix. `sample.num` realizations sampled from the final pool
  (columns subset of `series.sim`).

- filter_summary:

  Data frame summarizing pass counts and pass rates by criterion, and
  the final selection mode (tiered relaxation vs fallback).

- diagnostics:

  List containing window metadata, pool/sample indices, relaxation log,
  and final parameter values used for filtering.

- plots:

  NULL or a named list of ggplot objects when `make.plots = TRUE`.

## Details

**Length harmonization**: The evaluation is performed over
`n_use = min(length(series.obs), nrow(series.sim))`. If the observed
series is longer than `n_use`, one contiguous window of length `n_use`
is sampled and used for all comparisons. If simulated realizations are
longer than `n_use`, each realization is evaluated on its own sampled
contiguous window of length `n_use`. Returned series in `subsetted` and
`sampled` always use the original rows from `series.sim` (no trimming in
the outputs).

**Tail-mass metrics**: Lower and upper thresholds are defined by the
observed quantiles at `bounds$tail.low.p` and `bounds$tail.high.p`. Tail
behavior is summarized as normalized deficit/excess mass relative to
those thresholds, and compared using a stabilized log-distance with
`bounds$tail.eps`.

**Wavelet filter (observed-relevant GWS)**: The observed GWS is tested
against its significance curve to identify significant periods, which
are grouped into contiguous regions. For each realization, regional
integrated power and region-to-background contrast are compared to the
observed within the specified tolerances, requiring at least
`bounds$sig.frac` of regions to pass. Optional presence checks ensure
simulated power exceeds a minimum fraction of observed power in the
observed-relevant band.

## See also

[`wavelet_spectral_analysis`](https://github.com/Deltares-research/weathergenr/reference/wavelet_spectral_analysis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
out <- filter_warm_simulations(
  series.obs = obs,
  series.sim = sim_mat,
  sample.num = 5,
  seed = 123,
  make.plots = TRUE,
  verbose = TRUE
)
out$filter_summary
out$plots$wavelet_gws
} # }
```
