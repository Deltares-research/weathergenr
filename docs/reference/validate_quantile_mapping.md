# Validate precipitation quantile mapping adjustments

Convenience wrapper around \`diagnose_precip_qm()\` with argument names
aligned to QM workflows.

## Usage

``` r
validate_quantile_mapping(
  precip_org,
  precip_adjusted,
  month = NULL,
  year = NULL,
  mean_factor = NULL,
  var_factor = NULL,
  wet_thresh = 0.1,
  probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
)
```

## Arguments

- precip_org:

  Numeric vector. Original precipitation values.

- precip_adjusted:

  Numeric vector. Adjusted precipitation values (same length as
  \`precip_org\`).

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

A \`precip_qm_diagnostics\` object (see \`diagnose_precip_qm()\`).
