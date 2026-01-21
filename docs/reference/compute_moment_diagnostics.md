# Compute moment diagnostics (scenario-aware when factors + month/year are available)

Compute moment diagnostics (scenario-aware when factors + month/year are
available)

## Usage

``` r
compute_moment_diagnostics(
  precip_ref,
  precip_adj,
  mask,
  month = NULL,
  year = NULL,
  mean_factor = NULL,
  var_factor = NULL,
  tol_mean = 0.05,
  tol_var = 0.1
)
```
