# Compute Diagnostics for Precipitation Quantile Mapping

Calculates validation-style diagnostics comparing reference and
quantile-mapped precipitation series. Evaluates moment changes, quantile
preservation, tail behavior, dry/wet-day frequencies, spell lengths, and
optional temporal/monthly patterns.

## Usage

``` r
diagnose_prcp_qm(
  prcp_ref,
  prcp_adj,
  month = NULL,
  year = NULL,
  mean_factor = NULL,
  var_factor = NULL,
  wet_thresh = 0.1,
  probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
)
```

## Arguments

- prcp_ref:

  Numeric vector. Reference (original) precipitation values.

- prcp_adj:

  Numeric vector. Adjusted precipitation values (same length as
  `prcp_ref`).

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

A list of class `"prcp_qm_diagnostics"` containing:

- `moments`: Data frame of moment statistics

- `quantiles`: Data frame of quantile comparisons

- `extremes`: Data frame of extreme value metrics

- `temporal`: Data frame of temporal pattern metrics (optional)

- `monthly`: Data frame of monthly diagnostics (optional)

- `spells`: Data frame of spell length statistics

- `drydays`: Data frame of dry/wet-day frequency diagnostics

- `summary`: Overall assessment metrics
