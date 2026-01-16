# Standardize obs and sim to the same full-year length via random windowing

Applies the chosen window consistently across all grids and all
realizations after removing leap days. The window is drawn from the
longest contiguous full-year blocks available in each series.

## Usage

``` r
.align_obs_sim_periods(daily_obs, daily_sim, n_realizations, variables)
```

## Arguments

- daily_obs:

  List of observed data frames (one per grid).

- daily_sim:

  List of simulated realizations; each realization is a list of grids.

- n_realizations:

  Integer number of realizations.

- variables:

  Character vector of variables used for evaluation.

## Value

List with standardized \`daily_obs\`, \`daily_sim\`, and window
metadata.
