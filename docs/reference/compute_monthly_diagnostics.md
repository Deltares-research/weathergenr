# Compute month-by-month diagnostics (scenario-aware when factors + year are available)

Compute month-by-month diagnostics (scenario-aware when factors + year
are available)

## Usage

``` r
compute_monthly_diagnostics(
  precip_ref,
  precip_adj,
  month,
  mean_factor = NULL,
  var_factor = NULL,
  year = NULL,
  tol_mean = 0.05,
  tol_var = 0.1
)
```
