# Evaluate Stochastic Weather Generator Performance

Comprehensive diagnostic evaluation comparing synthetic weather
simulations against historical observations across multiple grid cells.
Computes summary statistics, wet/dry day counts, spell lengths, and
inter-site correlations, and generates diagnostic plots to assess
stochastic weather generator performance.

## Usage

``` r
evaluate_weather_generator(
  daily.sim = NULL,
  daily.obs = NULL,
  variables = NULL,
  variable.labels = NULL,
  realization.num = NULL,
  wet.quantile = 0.2,
  extreme.quantile = 0.8,
  output.path = NULL,
  save.plots = TRUE,
  show.title = TRUE,
  max.grids = 25,
  seed = NULL
)
```

## Arguments

- daily.sim:

  List of simulated weather realizations. Each element should be a list
  of data frames (one per grid cell), containing daily values and a
  \`date\` column.

- daily.obs:

  List of observed weather data frames (one per grid cell). Each should
  contain a \`date\` column and the variables specified in
  \`variables\`.

- variables:

  Character vector of variable names to evaluate (e.g., \`c("precip",
  "temp")\`). Must include "precip" for wet/dry spell analysis.

- variable.labels:

  Optional character vector of variable labels for plots. Defaults to
  \`variables\` if \`NULL\`.

- realization.num:

  Integer. Number of synthetic realizations in \`daily.sim\`.

- wet.quantile:

  Numeric between 0 and 1. Quantile threshold for wet days (default =
  0.2).

- extreme.quantile:

  Numeric between 0 and 1. Quantile threshold for extremely wet days
  (default = 0.8).

- output.path:

  Character. Directory path to save generated plots. If \`NULL\`, plots
  are not saved to disk.

- save.plots:

  Logical. Whether to save plots to \`output.path\` (default =
  \`TRUE\`).

- show.title:

  Logical. Whether to display titles in plots (default = \`TRUE\`).

- max.grids:

  Integer. Maximum number of grid cells to evaluate (default = 25). If
  more grids are provided, a random subsample is used to control memory.

- seed:

  Optional integer. Random seed for reproducible grid subsampling and
  year window selection. If NULL, results will vary between runs.

## Value

A named list of \`ggplot2\` plot objects with class
"weather_assessment". The returned object also contains fit summary
metrics stored as attributes.
