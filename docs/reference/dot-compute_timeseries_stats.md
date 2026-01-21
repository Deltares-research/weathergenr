# Compute time series statistics

Compute time series statistics

## Usage

``` r
.compute_timeseries_stats(data, variables, mc_thresholds)
```

## Arguments

- data:

  Data frame with daily values and a \`date\` column.

- variables:

  Character vector of variable names to summarize.

- mc_thresholds:

  Data frame of wet-day thresholds by grid and month.

## Value

List of seasonal, monthly, and annual stats plus wet/dry and correlation
data.
