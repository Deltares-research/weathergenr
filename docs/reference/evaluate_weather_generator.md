# Evaluate Stochastic Weather Generator Performance

Run a comprehensive diagnostic evaluation comparing synthetic weather
simulations against historical observations across multiple grid cells.
Computes summary statistics, wet/dry day counts, spell lengths, and
inter-site correlations, and generates diagnostic plots to assess
stochastic weather generator performance.

## Usage

``` r
evaluate_weather_generator(
  daily_sim = NULL,
  daily_obs = NULL,
  variables = NULL,
  variable_labels = NULL,
  n_realizations = NULL,
  wet_quantile = 0.2,
  extreme_quantile = 0.8,
  output_path = NULL,
  save_plots = TRUE,
  show_title = TRUE,
  verbose = TRUE,
  max_grids = 25,
  seed = NULL
)
```

## Arguments

- daily_sim:

  List of simulated weather realizations. Each element should be a list
  of data frames (one per grid cell), containing daily values and a
  \`date\` column.

- daily_obs:

  List of observed weather data frames (one per grid cell). Each should
  contain a \`date\` column and the variables specified in
  \`variables\`.

- variables:

  Character vector of variable names to evaluate (e.g., \`c("precip",
  "temp")\`). Must include "precip" for wet/dry spell analysis.

- variable_labels:

  Optional character vector of variable labels for plots. Defaults to
  \`variables\` if \`NULL\`.

- n_realizations:

  Integer. Number of synthetic realizations in \`daily_sim\`.

- wet_quantile:

  Numeric between 0 and 1. Quantile threshold for wet days (default =
  0.2).

- extreme_quantile:

  Numeric between 0 and 1. Quantile threshold for extremely wet days
  (default = 0.8).

- output_path:

  Character. Directory path to save generated plots. If \`NULL\`, plots
  are not saved to disk.

- save_plots:

  Logical. Whether to save plots to \`output_path\` (default =
  \`TRUE\`).

- show_title:

  Logical. Whether to display titles in plots (default = \`TRUE\`).

- verbose:

  Logical. Whether to emit console messages and the fit summary table.

- max_grids:

  Integer. Maximum number of grid cells to evaluate (default = 25). If
  more grids are provided, a random subsample is used to control memory.

- seed:

  Optional integer. Random seed for reproducible grid subsampling and
  year window selection. If NULL, results will vary between runs.

## Value

A named list of \`ggplot2\` plot objects with class
"weather_assessment". The returned object also contains attributes:

- `fit_summary`: data frame of per-realization fit metrics and ranks.

- `metadata`: list with `n_grids`, `n_realizations`, `variables`, and
  `assessment_date`.

## Details

The function standardizes simulated and observed series to full-year
windows (after leap-day removal), optionally subsamples grid cells to
manage memory, and returns diagnostic plots alongside a summarized fit
table. Missing values are ignored in summary statistics (\`na.rm =
TRUE\`), and correlations use pairwise complete observations. Use
\`seed\` for reproducible subsampling and window selection; the original
RNG state is restored on exit.

## Examples

``` r
set.seed(1)
dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
dates <- dates[format(dates, "%m-%d") != "02-29"]
obs_grid <- list(data.frame(
  date = dates,
  precip = rgamma(length(dates), shape = 2, scale = 2),
  temp = rnorm(length(dates), mean = 10, sd = 3)
))
sim_grid <- list(list(data.frame(
  date = dates,
  precip = rgamma(length(dates), shape = 2, scale = 2),
  temp = rnorm(length(dates), mean = 10, sd = 3)
)))
out <- evaluate_weather_generator(
  daily_sim = sim_grid,
  daily_obs = obs_grid,
  variables = c("precip", "temp"),
  n_realizations = 1,
  output_path = NULL,
  save_plots = FALSE,
  show_title = FALSE
)
#> INFO [2026-01-16 18:36:41] [VALIDATE] Start | grids = 1 | realizations = 1
#> INFO [2026-01-16 18:36:41] [VALIDATE] Variables = precip,temp
#> INFO [2026-01-16 18:36:41] [VALIDATE] Parameters: wet.q = 0.2 | extreme.q = 0.8
#> INFO [2026-01-16 18:36:41] [VALIDATE] Standardizing obs/sim periods to full years and equal length
#> INFO [2026-01-16 18:36:41] [VALIDATE] Standardized period | Obs = 2001-2001 | Sim = 2001-2001
#> INFO [2026-01-16 18:36:41] [VALIDATE] Processing observed data
#> INFO [2026-01-16 18:36:41] [VALIDATE] Processing simulated data (1 realizations)
#> INFO [2026-01-16 18:36:41] [VALIDATE] Preparing diagnostic data for plotting
#> INFO [2026-01-16 18:36:41] [VALIDATE] Generating diagnostic plots
#> Warning: There were 2 warnings in `dplyr::summarise()`.
#> The first warning was:
#> ℹ In argument: `.min = min(c(.data[["Observed"]], .data[["Simulated"]]), na.rm
#>   = TRUE)`.
#> Caused by warning in `min()`:
#> ! no non-missing arguments to min; returning Inf
#> ℹ Run dplyr::last_dplyr_warnings() to see the 1 remaining warning.
#> INFO [2026-01-16 18:36:42] [VALIDATE] Generated 11 diagnostic plots.
#> INFO [2026-01-16 18:36:42] [VALIDATE] Computing fit metrics for all realizations
#> Warning: There were 2 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `dplyr::across(...)`.
#> Caused by warning in `min()`:
#> ! no non-missing arguments to min; returning Inf
#> ℹ Run dplyr::last_dplyr_warnings() to see the 1 remaining warning.
#> INFO [2026-01-16 18:36:42] [VALIDATE] Displaying fit assessment summary
#> 
#> =============================================================================================== 
#>  FIT ASSESSMENT SUMMARY - ALL REALIZATIONS
#> =============================================================================================== 
#> 
#>  Metrics shown (Mean Absolute Error, lower is better):
#>   MAE measures average magnitude of errors between simulated and observed
#> 
#>   - Mean.precip   : MAE of monthly means
#>   - SD.precip     : MAE of monthly standard deviations
#>   - Days.Wet   : MAE of wet/dry day counts
#>   - Spell.Wet  : MAE of wet/dry spell lengths (days)
#>   - Cor.Cross   : MAE of cross-grid correlations
#>   - Cor.Inter   : MAE of inter-variable correlations
#>   - Score       : Normalized overall score across all metrics
#> 
#>    Rlz   Rank   Mean.precip   SD.precip   Days.Wet   Spell.Wet   Cor.Cross   Cor.Inter    Score 
#> ----------------------------------------------------------------------------------------------- 
#>      1      1        1.6773      1.1552    11.5091     10.2257          NA      0.0479   0.0000 
#> =============================================================================================== 
#> 
#>  Summary:
#>   - Best realization  : 1 (score = 0.0000)
#>   - Worst realization : 1 (score = 0.0000)
#>   - Median score      : 0.0000
#> 
#> INFO [2026-01-16 18:36:42] [VALIDATE] Assessment completed successfully
class(out)
#> [1] "weather_assessment" "list"              
```
