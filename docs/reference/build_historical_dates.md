# Build Historical Date Table

Constructs a tibble of historical dates with calendar and water-year
indices, then restricts to complete years (365 days each) and returns
the longest contiguous block of complete years.

## Usage

``` r
build_historical_dates(obs_dates, year_start_month, verbose = FALSE)
```

## Arguments

- obs_dates:

  Date vector of observation dates (leap days already removed).

- year_start_month:

  Integer 1-12. First month of the simulation year.

- verbose:

  Logical. If TRUE, logs summary information.

## Value

A list with three elements:

- dates_df:

  Tibble with columns: dateo, year, wyear, month, day, date.

- wyear_idx:

  Integer vector of row indices in obs_dates corresponding to complete
  years.

- complete_wyears:

  Integer vector of complete water years.
