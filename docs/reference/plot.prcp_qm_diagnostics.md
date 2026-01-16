# Plot diagnostics for precipitation QM

Plots selected diagnostic panels from a `"prcp_qm_diagnostics"` object.

## Usage

``` r
# S3 method for class 'prcp_qm_diagnostics'
plot(x, which = c("all", "moments", "quantiles", "extremes", "monthly"), ...)
```

## Arguments

- x:

  A `"prcp_qm_diagnostics"` object.

- which:

  Character. Which plot to return.

- ...:

  Additional arguments passed to plotting functions (unused).

## Value

A ggplot object when a single plot is requested, or a list of plots when
`which = "all"` and `gridExtra` is unavailable.
