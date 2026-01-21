# Compute Area-Averaged Daily and Annual Climate Series

Computes area-averaged (mean across grid cells) daily climate values and
aggregates them to annual means by water year.

## Usage

``` r
compute_area_averages(obs_data, wyear_idx, wyear, vars)
```

## Arguments

- obs_data:

  Named list of data frames (one per grid cell).

- wyear_idx:

  Integer vector of row indices to extract from each data frame.

- wyear:

  Integer vector of water years corresponding to wyear_idx.

- vars:

  Character vector of variable names to average.

## Value

A list with two elements:

- daily:

  Data frame of area-averaged daily values with columns for each
  variable plus wyear.

- annual:

  Tibble of annual means with columns wyear plus each variable.
