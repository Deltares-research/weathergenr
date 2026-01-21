# Prepare Generator Output for Evaluation

Transforms the output from
[`generate_weather`](https://deltares-research.github.io/weathergenr/reference/generate_weather.md)
into the format required by
[`evaluate_weather_generator`](https://deltares-research.github.io/weathergenr/reference/evaluate_weather_generator.md).
This helper function handles the date resampling, complete-year
filtering, and data extraction that is typically needed between
generation and evaluation steps.

## Usage

``` r
prepare_evaluation_data(
  gen_output,
  obs_data,
  obs_dates,
  grid_ids,
  variables,
  min_days_per_year = 365L,
  verbose = TRUE
)
```

## Arguments

- gen_output:

  List returned by
  [`generate_weather`](https://deltares-research.github.io/weathergenr/reference/generate_weather.md),
  containing `resampled` (tibble of resampled observation dates per
  realization) and `dates` (vector of simulated dates).

- obs_data:

  Named list of data frames from
  [`read_netcdf`](https://deltares-research.github.io/weathergenr/reference/read_netcdf.md)`$data`,
  one per grid cell. Each data frame must contain columns for all
  variables.

- obs_dates:

  Date vector corresponding to rows in each `obs_data` element.
  Typically from `read_netcdf()$date`.

- grid_ids:

  Character or integer vector of grid cell identifiers to include in the
  evaluation. Must match names or indices in `obs_data`.

- variables:

  Character vector of variable names to extract (e.g.,
  `c("precip", "temp")`).

- min_days_per_year:

  Integer. Minimum number of days required to consider a year complete.
  Default is 365. Use 360 for 360-day calendars.

- verbose:

  Logical. If `TRUE`, prints progress messages to console via the
  internal
  [`.log()`](https://deltares-research.github.io/weathergenr/reference/dot-log.md)
  function. Default is `TRUE`.

## Value

A list with two elements:

- sim_data:

  List of length `n_realizations`. Each element is a list of data frames
  (one per grid cell) with columns `date` followed by the requested
  `variables`.

- obs_data:

  List of data frames (one per grid cell) with columns `date` followed
  by the requested `variables`.

## Details

The function performs the following transformations:

1.  Maps resampled dates back to row indices in the original
    observations

2.  Filters the simulation period to complete years only (\>= 365 days)

3.  Extracts the specified variables from observed data using resampled
    indices

4.  Formats both simulated and observed data as required by the
    evaluator

Complete years are identified by calendar year boundaries. For
water-year simulations, users should ensure the simulation spans full
water years.

## See also

[`generate_weather`](https://deltares-research.github.io/weathergenr/reference/generate_weather.md),
[`evaluate_weather_generator`](https://deltares-research.github.io/weathergenr/reference/evaluate_weather_generator.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running generate_weather()
ncdata <- read_netcdf("climate_data.nc")
gen_output <- generate_weather(
  obs_data = ncdata$data,
  obs_grid = ncdata$grid,
  obs_dates = ncdata$date,
  ...
)

# Prepare data for evaluation
eval_data <- prepare_evaluation_data(
  gen_output = gen_output,
  obs_data   = ncdata$data,
  obs_dates  = ncdata$date,
  grid_ids   = ncdata$grid$id,
  variables  = c("precip", "temp")
)

# Run evaluation
results <- evaluate_weather_generator(
  daily_sim = eval_data$sim_data,
  daily_obs = eval_data$obs_data,
  ...
)
} # }
```
