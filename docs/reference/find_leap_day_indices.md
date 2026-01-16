# Find leap-day (Feb 29) indices in a date vector

Scans a vector of dates and returns the integer indices of all elements
that correspond to February 29. If the vector contains no leap days,
returns `NULL`.

## Usage

``` r
find_leap_day_indices(date)
```

## Arguments

- date:

  A vector coercible to `Date`.

## Value

Integer vector of indices where dates equal February 29, or `NULL` if
none.

## Examples

``` r
find_leap_day_indices(as.Date(c("1980-02-28", "1980-02-29", "1981-01-01")))
#> [1] 2
find_leap_day_indices(as.Date(c("2001-01-01", "2001-12-31")))
#> NULL
```
