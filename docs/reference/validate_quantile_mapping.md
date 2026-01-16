# Validate precipitation quantile mapping adjustments

Computes validation diagnostics comparing original and adjusted
precipitation series using the same metrics returned by
[`diagnose_prcp_qm`](https://deltares-research.github.io/weathergenr/reference/diagnose_prcp_qm.md).

## Usage

``` r
validate_quantile_mapping(
  prcp_org,
  prcp_adjusted,
  month = NULL,
  year = NULL,
  mean_factor = NULL,
  var_factor = NULL,
  wet_thresh = 0.1,
  probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
)
```

## Arguments

- prcp_org:

  Numeric vector. Original precipitation values.

- prcp_adjusted:

  Numeric vector. Adjusted precipitation values (same length as
  `prcp_org`).

- month:

  Integer vector (optional). Month index (1–12) for each observation.

- year:

  Integer vector (optional). Year index (1..n_years) for each
  observation.

- mean_factor:

  Numeric matrix (optional). Target mean scaling factors (n_years x 12).

- var_factor:

  Numeric matrix (optional). Target variance scaling factors (n_years x
  12).

- wet_thresh:

  Numeric. Threshold for defining wet days (default = 0.1 mm).

- probs:

  Numeric vector. Quantile probabilities to evaluate.

## Value

A list of class `"prcp_qm_diagnostics"` with diagnostic tables and
summary information (see
[`diagnose_prcp_qm`](https://deltares-research.github.io/weathergenr/reference/diagnose_prcp_qm.md)
for details).

## Details

This is a convenience wrapper around
[`diagnose_prcp_qm`](https://deltares-research.github.io/weathergenr/reference/diagnose_prcp_qm.md)
with argument names aligned to quantile-mapping workflows.

## Examples

``` r
prcp_org <- c(0, 1, 2, 0, 5, 3)
prcp_adjusted <- c(0, 1.1, 2.2, 0, 4.8, 3.1)
month <- c(1, 1, 1, 1, 1, 1)
year <- c(1, 1, 1, 1, 1, 1)
validate_quantile_mapping(prcp_org, prcp_adjusted, month = month, year = year)
#> 
#> === Precipitation QM Diagnostics ===
#> 
#> Overall Assessment: excellent 
#> Overall Score: 83.2 / 100
#> 
#> Component Scores:
#>   moments: 20.0
#>   quantiles: 24.0
#>   extremes: 19.2
#>   spells: 10.0
#>   drydays: 10.0
#> 
#> Key Findings:
#>   Mean change: 1.8%
#>   Variance change: -16.1%
#>   Dry days preserved: TRUE
#>   Extreme amplification: 0.96x
#> 
#> Moment Preservation:
#>    metric original adjusted pct_change assessment
#>      mean    2.750    2.800      1.818  excellent
#>        sd    1.708    1.564     -8.411       good
#>  variance    2.917    2.447    -16.114       poor
#>        cv    0.621    0.559    -10.046       good
#>  skewness    0.282    0.189    -32.943       poor
#>  kurtosis   -1.962   -1.977      0.788       good
#> 
#> Dry/Wet Day Preservation:
#>  category original adjusted assessment
#>  dry_days        2        2    perfect
#>  wet_days        4        4    perfect
#> 
```
