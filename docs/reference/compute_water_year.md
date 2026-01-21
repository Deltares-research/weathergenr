# Calculate water year from a date vector

Assigns a water year to each date based on a custom starting month. A
water year groups months across calendar years (e.g., Oct-Sep).

## Usage

``` r
compute_water_year(date, water_year_start_month = 1)
```

## Arguments

- date:

  A vector of `Date` objects.

- water_year_start_month:

  Integer (1-12) indicating the first month of the water year. For
  example, `10` for October or `6` for June. Default is `1` (calendar
  year).

## Value

An integer vector of the same length as `date`, giving the water year of
each date.

## Examples

``` r
dates <- as.Date(c("2022-09-30", "2022-10-01", "2023-06-15"))
compute_water_year(dates, water_year_start_month = 10)
#> [1] 2022 2023 2023
```
