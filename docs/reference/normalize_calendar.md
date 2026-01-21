# Normalize Calendar to 365 Days

Removes February 29 (leap days) from observation dates and corresponding
rows from all grid cell data frames to enforce a 365-day calendar.

## Usage

``` r
normalize_calendar(obs_dates, obs_data, verbose = FALSE)
```

## Arguments

- obs_dates:

  Date vector of observation dates.

- obs_data:

  Named list of data frames (one per grid cell).

- verbose:

  Logical. If TRUE, logs the number of dropped days.

## Value

A list with two elements:

- dates:

  Date vector with leap days removed.

- data:

  List of data frames with corresponding rows removed.
