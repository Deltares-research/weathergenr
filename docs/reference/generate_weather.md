# Generate Synthetic Gridded Daily Weather Series

Generates stochastic gridded daily weather by coupling an annual-scale
low-frequency generator with daily-scale resampling and persistence
logic. The workflow combines:

- Wavelet Autoregressive Modeling (WARM) on an annual aggregate of
  `warm_var` to simulate low-frequency variability,

- annual K-nearest-neighbor (KNN) matching to select historical analogue
  years,

- a three-state (dry, wet, extreme) daily Markov chain to control spell
  persistence,

- daily KNN resampling of precipitation and temperature anomalies.

The simulation-year definition is inferred from `year_start_month`:

- `year_start_month == 1`: calendar-year simulation,

- `year_start_month != 1`: water-year simulation starting in
  `year_start_month`.

## Usage

``` r
generate_weather(
  obs_data = NULL,
  obs_grid = NULL,
  obs_dates = NULL,
  vars = NULL,
  n_years = NULL,
  start_year = 2020,
  year_start_month = 1,
  n_realizations = 5,
  warm_var = "precip",
  warm_signif = 0.9,
  warm_pool_size = 5000,
  warm_filter_bounds = list(),
  relax_priority = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
  annual_knn_n = 120,
  wet_q = 0.3,
  extreme_q = 0.8,
  dry_spell_factor = rep(1, 12),
  wet_spell_factor = rep(1, 12),
  out_dir = tempdir(),
  seed = NULL,
  parallel = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  save_plots = TRUE
)
```

## Arguments

- obs_data:

  Named list of data.frames (one per grid cell) with observed daily
  weather. Each data.frame must have one row per day corresponding to
  `obs_dates`, and columns for all `vars`.

- obs_grid:

  Data.frame describing grid cells. Must contain columns `xind`, `yind`,
  `x`, `y`. If `id` is missing, it is created as a sequential index.

- obs_dates:

  Vector of class `Date` corresponding to rows of each `obs_data[[i]]`.
  May include Feb 29; Feb 29 is removed internally along with the
  corresponding rows in `obs_data`.

- vars:

  Character vector of daily variables to simulate (e.g.,
  `c("precip","temp")`).

- n_years:

  Integer number of years to simulate. If `NULL`, defaults to the number
  of complete historical simulation-years available after enforcing the
  365-day calendar and restricting to the longest contiguous block of
  complete years.

- start_year:

  Integer first simulation year (calendar year if
  `year_start_month == 1`, otherwise first water year label).

- year_start_month:

  Integer in `1:12`. First month of the simulation year. Use `1` for
  calendar years; use `10` for Oct–Sep water years, etc.

- n_realizations:

  Integer number of realizations to generate.

- warm_var:

  Character name of the annual driver variable used in WARM. Must be
  present in `vars`.

- warm_signif:

  Numeric in (0,1). Significance level for the wavelet spectrum tests,
  against a red-noise background.

  It does **not** select which wavelet components are simulated: every
  MODWT MRA component is modelled, significant or not. See the WARM
  section of *Details*. Its two live effects are:

  - the threshold above which an observed spectral peak counts as
    significant in
    [`filter_warm_pool`](https://deltares-research.github.io/weathergenr/reference/filter_warm_pool.md),
    and therefore must be reproduced by a candidate trace to pass the
    `wavelet` criterion;

  - raising the MODWT level count when the continuous wavelet transform
    finds a significant period. On records of 20-40 years this test
    rarely fires, so in practice the level count is set by series length
    alone.

- warm_pool_size:

  Integer number of candidate annual traces to generate before filtering
  down to `n_realizations`.

- warm_filter_bounds:

  Named list of filtering thresholds and relaxation controls forwarded
  to
  [`filter_warm_pool`](https://deltares-research.github.io/weathergenr/reference/filter_warm_pool.md)
  as `filter_bounds`. Any entry overrides internal defaults. Uses
  snake_case keys (e.g. `tail_low_p`).

- relax_priority:

  Character vector giving the relaxation priority for WARM filtering,
  forwarded to
  [`filter_warm_pool`](https://deltares-research.github.io/weathergenr/reference/filter_warm_pool.md)
  as `relax_order`. Must contain each of
  `c("mean","sd","tail_low","tail_high","wavelet")` exactly once.
  Filters are relaxed iteratively by loosening the currently most
  restrictive criterion (lowest pass rate), subject to this priority
  ordering.

- annual_knn_n:

  Integer number of historical analogue years considered in the annual
  KNN step of the disaggregation.

- wet_q:

  Numeric in (0,1). Wet-day threshold quantile used in the daily Markov
  chain classification.

- extreme_q:

  Numeric in (0,1). Extreme-day threshold quantile used in the daily
  Markov chain classification. Must be greater than `wet_q`.

- dry_spell_factor:

  Numeric length-12 vector. Per-month multiplier on the **mean dry spell
  length**: 1.5 makes dry spells half again as long, 2 doubles them. 1
  (the default) leaves the month untouched.

- wet_spell_factor:

  Numeric length-12 vector. The same, for mean wet spell length, where a
  wet spell is an unbroken run of wet or extreme days.

  Note these are not independent of wet-day frequency: in a Markov
  chain, time spent dry and time spent wet must trade off, so
  lengthening dry spells alone lowers the wet-day fraction (on the
  packaged example, a `dry_spell_factor` of 2 lowers it by about 10
  percent). That is a property of the chain rather than an artefact.
  Adjust both factors together to lengthen spells while holding the
  wet-day fraction roughly fixed.

- out_dir:

  Character directory for diagnostics and CSV outputs. The directory is
  created if it does not exist.

- seed:

  Optional integer seed for reproducibility. Used to initialize RNG
  streams for annual (WARM) and daily resampling components.

- parallel:

  Logical; if `TRUE`, run the daily disaggregation across realizations
  in parallel using a PSOCK cluster.

- n_cores:

  Optional integer number of cores for parallel execution. If `NULL` and
  `parallel = TRUE`, defaults to `max(1, parallel::detectCores() - 1)`.

- verbose:

  Logical; if `TRUE`, emits progress logs via
  [`.log()`](https://deltares-research.github.io/weathergenr/reference/dot-log.md).

- save_plots:

  Logical; if `TRUE`, writes plot to `output_dir`.

## Value

A list with:

- resampled:

  Tibble/data.frame with one column per realization. Each column
  contains historical observation dates (`dateo`) selected as analogues
  for each simulated day. Column names are `rlz_1`, `rlz_2`, ...

- dates:

  Vector of class `Date` giving the simulated daily time axis (Feb 29
  excluded).

## Details

**Calendar handling (robust Gregorian input):** Inputs may be on the
Gregorian calendar and may include Feb 29. Internally, the function
enforces a 365-day calendar by removing Feb 29 from `obs_dates` and
dropping the corresponding rows from every element of `obs_data`. All
downstream processing and outputs therefore use a 365-day calendar.

**Workflow:**

1.  **Historical preprocessing:** Enforce a 365-day calendar, compute
    calendar/water-year indices, and build a historical dates table used
    for KNN matching and Markov chain sequencing.

2.  **Annual WARM simulation:** Aggregate historical daily data to
    annual means (by water year if applicable), run wavelet analysis on
    the annual series, simulate candidate annual traces via
    [`simulate_warm`](https://deltares-research.github.io/weathergenr/reference/simulate_warm.md),
    and subset traces using
    [`filter_warm_pool`](https://deltares-research.github.io/weathergenr/reference/filter_warm_pool.md).

3.  **Daily disaggregation:** For each realization, call
    [`resample_weather_dates`](https://deltares-research.github.io/weathergenr/reference/resample_weather_dates.md)
    to generate a simulated sequence of historical analogue dates (in
    the internal 365-day calendar).

4.  **Output construction:** Map resampled internal dates back to
    historical observation dates (`dateo`) and return simulated dates.

**WARM, and how it relates to Steinschneider & Brown (2013):** The
annual generator follows the wavelet autoregressive model (WARM) of Kwon
et al. (2007), as adopted by Steinschneider & Brown (2013), but differs
in the decomposition and in component selection. Both differences are
deliberate; neither is a partial implementation.

In the original, a continuous wavelet transform (CWT) decomposes the
annual series, the scales whose global wavelet power exceeds a red-noise
significance level are retained as low-frequency signals, and the model
is an autoregression per retained signal *plus* an autoregression on the
residual. Nothing is discarded: what is not retained becomes the
residual term. Steinschneider & Brown's own application retained a
single signal.

This package decomposes with a MODWT multiresolution analysis instead.
That basis is exactly additive – the detail bands and the smooth sum
back to the input, with reconstruction error at machine precision – so
there is no residual left over to model separately, and every component
is fitted with an ARMA. On a 20-year annual series the decomposition
yields two components, `D1` and `S1`, so the fitted structure reduces to
the same shape as the original with one retained signal. That is a
consequence of how few levels a short record supports, not a
correspondence between a MODWT band and the residual term: a longer
record yields `D1`, `D2`, `S2` and no band plays that role.

Significance is therefore *not* used to select components, for two
reasons. There is no residual term to relegate a non-significant band
to, so dropping one would discard variance rather than reclassify it –
on the packaged example `D1` carries 23 percent of the total. And on
records of this length the CWT global-power test against a red-noise
background detects nothing at any conventional level, so a selection
rule would either never fire or would collapse the model to a single
ARMA on the whole series, discarding the scale structure the method
exists to represent.

Note that MODWT MRA components are additive but not orthogonal, so their
variances do not sum to the total; on the packaged example the
cross-covariance accounts for 18 percent.
[`simulate_warm`](https://deltares-research.github.io/weathergenr/reference/simulate_warm.md)
corrects for this when recombining simulated components.

Daily disaggregation currently requires `vars` to include `"precip"` and
`"temp"`, which are used as precipitation and temperature drivers in
KNN + Markov resampling.

## References

Kwon, H.-H., Lall, U., & Khalil, A. F. (2007). Stochastic simulation
model for nonstationary time series using an autoregressive wavelet
decomposition: Applications to rainfall and temperature. *Water
Resources Research*, 43(5), W05407.
[doi:10.1029/2006WR005258](https://doi.org/10.1029/2006WR005258)

Steinschneider, S., & Brown, C. (2013). A semiparametric multivariate,
multisite weather generator with low-frequency variability for use in
climate risk assessments. *Water Resources Research*, 49(11), 7205-7220.
[doi:10.1002/wrcr.20528](https://doi.org/10.1002/wrcr.20528)

## Examples

``` r
if (FALSE) { # \dontrun{
# obs_data: list of grid-cell data.frames with columns precip/temp
# obs_grid: data.frame with xind/yind/x/y
# obs_dates: Date vector aligned with rows of each obs_data[[i]]
out <- generate_weather(
  obs_data = obs_data,
  obs_grid = obs_grid,
  obs_dates = obs_dates,
  vars = c("precip", "temp"),
  start_year = 2020,
  n_years = 30,
  year_start_month = 10,
  n_realizations = 20,
  warm_var = "precip",
  warm_signif = 0.90,
  warm_pool_size = 5000,
  annual_knn_n = 120,
  wet_q = 0.30,
  extreme_q = 0.80,
  out_dir = tempdir(),
  seed = 123,
  parallel = TRUE,
  n_cores = 4,
  verbose = TRUE
)
} # }
```
