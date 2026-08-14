# Row index selecting the rows of a shared date vector falling in given years

Every grid within \`daily_obs\` (and within one realization of
\`daily_sim\`) shares a single date vector, so the row index is computed
once here and reused across grids rather than recomputing the year
vector per grid. With 25 grids and 3 realizations that is one call
instead of 100.

## Usage

``` r
.year_keep_index(dates, years_keep, year_start_month = 1L)
```

## Arguments

- dates:

  Date vector shared by every grid on this side.

- years_keep:

  Integer vector of years (calendar or water years) to retain.

- year_start_month:

  Integer 1-12. When \> 1, uses water-year grouping.

## Value

Logical vector, length \`length(dates)\`.
