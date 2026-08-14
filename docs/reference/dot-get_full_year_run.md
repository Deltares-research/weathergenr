# Find full (365-day) years after leap-day removal and return longest contiguous block

Find full (365-day) years after leap-day removal and return longest
contiguous block

## Usage

``` r
.get_full_year_run(dates, year_start_month = 1L)
```

## Arguments

- dates:

  Date vector (leap days should already be removed).

- year_start_month:

  Integer 1-12. First month of the simulation year. When \> 1, groups by
  water year via \`compute_water_year()\`.

## Value

List with a \`years\` integer vector. Empty if no full-year blocks
exist.
