# Compute Diagnostics for Precipitation Quantile Mapping

Computes diagnostics comparing reference and quantile-mapped
precipitation series. Reports moment changes, quantiles, tail metrics,
dry/wet-day frequencies, spell lengths, and optional temporal/monthly
patterns. When \`mean_factor\`/\`var_factor\` are provided together with
\`month\` and \`year\`, moment and monthly diagnostics are compared
against the intended perturbations (scenario-aware null hypothesis).

## Usage

``` r
diagnose_precip_qm(
  precip_ref,
  precip_adj,
  month = NULL,
  year = NULL,
  mean_factor = NULL,
  var_factor = NULL,
  wet_thresh = 0.1,
  probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
)
```

## Arguments

- precip_ref:

  Numeric vector. Reference (original) precipitation values.

- precip_adj:

  Numeric vector. Adjusted precipitation values (same length as
  \`precip_ref\`).

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

A list of class \`precip_qm_diagnostics\` with elements:

- \`moments\`: moment diagnostics table

- \`quantiles\`: quantile comparison table

- \`extremes\`: tail metrics table

- \`temporal\`: temporal metrics table (optional)

- \`monthly\`: month-by-month diagnostics table (optional)

- \`spells\`: wet/dry spell diagnostics table

- \`drydays\`: wet/dry frequency diagnostics table

- \`summary\`: compact summary list
