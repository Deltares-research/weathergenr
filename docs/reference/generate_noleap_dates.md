# Generate a no-leap (365-day) daily date sequence

Constructs a continuous sequence of dates starting from a given start
date, using a no-leap calendar (i.e., February 29 removed). The function
repeats a 365-day month-day template and increments the year every 365
days.

## Usage

``` r
generate_noleap_dates(start_date, n_days)
```

## Arguments

- start_date:

  Date or character coercible to Date. The first date of the sequence.
  If the date is February 29, it is mapped to February 28 because the
  no-leap template does not include Feb 29.

- n_days:

  Integer. Total number of dates to generate.

## Value

A vector of class `Date` with length `n_days`, following a no-leap
(365-day) calendar.

## Examples

``` r
generate_noleap_dates("1980-01-01", 365)
generate_noleap_dates(as.Date("1995-02-28"), 1000)
```
