# Conditional precip-X correlations within each grid

Conditional precip-X correlations within each grid

## Usage

``` r
compute_conditional_precip_cor(
  data,
  variables,
  mc.thresholds = NULL,
  wet_def = c("gt0", "monthly_quantile"),
  use_anom = TRUE,
  use_log_precip_on_wet = TRUE,
  method = "pearson",
  min_pairs = 50
)
```
