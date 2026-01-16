# Generate Dummy Points to Enforce Symmetric Facet Axes

Creates invisible points with x = y at the min/max range per facet to
force ggplot to use symmetric x/y limits within each facet when scales
are free.

## Usage

``` r
generate_symmetric_dummy_points(df, facet_var, x_col, y_col)
```

## Arguments

- df:

  Data frame containing the plotted data.

- facet_var:

  Character; name of the facet column in `df`.

- x_col:

  Character; name of the x column in `df`.

- y_col:

  Character; name of the y column in `df`.

## Value

Data frame with columns `facet_var`, `Observed`, `Simulated` (using
`x_col` and `y_col` names) containing min/max dummy points.

## Details

The output includes two points per facet (min and max). These points can
be added with `geom_blank()` to enforce consistent axis limits without
affecting the visible data.

## Examples

``` r
df <- data.frame(
  variable = c("precip", "precip", "temp", "temp"),
  Observed = c(1, 5, 10, 12),
  Simulated = c(0.5, 6, 9, 13)
)
generate_symmetric_dummy_points(df, "variable", "Observed", "Simulated")
#> # A tibble: 4 × 3
#>   variable Observed Simulated
#>   <chr>       <dbl>     <dbl>
#> 1 precip        0.5       0.5
#> 2 temp          9         9  
#> 3 precip        6         6  
#> 4 temp         13        13  
```
