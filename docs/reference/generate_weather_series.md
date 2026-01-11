# Generate Synthetic Gridded Daily Weather Series

Generate stochastic gridded daily weather by coupling an annual-scale
low-frequency generator with daily-scale resampling and persistence
logic. The workflow combines:

- Wavelet Autoregressive Modeling (WARM) on an annual aggregate of
  `warm.variable` to simulate low-frequency variability,

- annual K-nearest-neighbor (KNN) matching to select historical analogue
  years,

- a three-state (dry, wet, extreme) daily Markov chain to control spell
  persistence,

- daily KNN resampling of precipitation and temperature anomalies.

The simulation regime is inferred from `month.start`:

- `month.start == 1`: calendar-year simulation,

- `month.start != 1`: water-year simulation starting at `month.start`.

## Usage

``` r
generate_weather_series(
  weather.data = NULL,
  weather.grid = NULL,
  weather.date = NULL,
  variables = NULL,
  sim.year.num = NULL,
  sim.year.start = 2020,
  month.start = 1,
  realization.num = 5,
  warm.variable = "precip",
  warm.signif.level = 0.9,
  warm.sample.num = 5000,
  knn.sample.num = 120,
  mc.wet.quantile = 0.3,
  mc.extreme.quantile = 0.8,
  dry.spell.change = rep(1, 12),
  wet.spell.change = rep(1, 12),
  output.path = tempdir(),
  seed = NULL,
  compute.parallel = FALSE,
  num.cores = NULL
)
```

## Arguments

- weather.data:

  Named list of data frames (one per grid cell) with observed daily
  weather. Each data frame must have one row per day corresponding to
  `weather.date`, and columns for all `variables`.

- weather.grid:

  Data frame describing grid cells. Must contain columns `id`, `xind`,
  `yind`, `x`, `y`.

- weather.date:

  Vector of class `Date` corresponding to rows of each
  `weather.data[[i]]`. May include Feb 29; Feb 29 is removed internally
  along with the corresponding rows in `weather.data`.

- variables:

  Character vector of daily variables to simulate (e.g.,
  `c("precip","temp")`).

- sim.year.num:

  Integer number of years to simulate. If `NULL`, defaults to the number
  of unique historical simulation-years (calendar year or water year
  depending on `month.start`).

- sim.year.start:

  Integer first simulation year (calendar year if `month.start == 1`,
  otherwise first water year).

- month.start:

  Integer in `1:12`. First month of the simulation year.

- realization.num:

  Integer number of realizations to generate.

- warm.variable:

  Character name of the annual driver variable used in WARM (must be in
  `variables`).

- warm.signif.level:

  Numeric in (0,1). Wavelet significance level for retaining
  low-frequency components in WARM.

- warm.sample.num:

  Integer number of candidate annual traces to generate before
  subsetting.

- knn.sample.num:

  Integer number of historical years sampled in annual KNN.

- mc.wet.quantile:

  Numeric in (0,1). Wet-day threshold quantile.

- mc.extreme.quantile:

  Numeric in (0,1). Extreme-day threshold quantile.

- dry.spell.change:

  Numeric length-12 vector of monthly dry-spell adjustment factors.

- wet.spell.change:

  Numeric length-12 vector of monthly wet-spell adjustment factors.

- output.path:

  Character directory for diagnostics/CSV outputs.

- seed:

  Optional integer seed for reproducibility.

- compute.parallel:

  Logical; if `TRUE`, run disaggregation in parallel.

- num.cores:

  Optional integer number of cores for parallel execution.

## Value

A list with:

- resampled:

  Tibble/data.frame with one column per realization. Each column
  contains the historical `dateo` values selected as analogues for each
  simulated day.

- dates:

  Vector of `Date` giving the simulated daily time axis (Feb 29
  excluded).

## Details

**Calendar handling (robust Gregorian input):** Inputs may be on the
Gregorian calendar and may include Feb 29. Internally, the function
enforces a 365-day calendar by removing Feb 29 from `weather.date` and
dropping the corresponding rows from every element of `weather.data`.
All downstream processing and outputs therefore use a 365-day calendar.

**Workflow:**

1.  **Historical preprocessing:** enforce a 365-day calendar, compute
    calendar/water-year indices, and build a historical dates table used
    for KNN matching and Markov chain sequencing.

2.  **Annual WARM simulation:** aggregate historical daily data to
    annual means (by water year if applicable), run wavelet analysis on
    the annual series, simulate candidate annual traces via
    [`wavelet_arima`](https://github.com/Deltares-research/weathergenr/reference/wavelet_arima.md),
    and subset traces using
    [`filter_warm_simulations`](https://github.com/Deltares-research/weathergenr/reference/filter_warm_simulations.md).

3.  **Daily disaggregation:** for each realization, call
    [`resample_weather_dates()`](https://github.com/Deltares-research/weathergenr/reference/resample_weather_dates.md)
    to generate a simulated sequence of historical analogue dates (in
    the internal 365-day calendar).

4.  **Output construction:** map resampled internal dates back to
    historical observation dates (`dateo`) and return simulated dates.
