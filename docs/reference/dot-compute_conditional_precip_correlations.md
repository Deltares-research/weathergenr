# Conditional precip-X correlations within each grid

Conditional precip-X correlations within each grid

## Usage

``` r
.compute_conditional_precip_correlations(
  data,
  variables,
  mc_thresholds = NULL,
  wet_def = c("gt0", "monthly_quantile"),
  use_anom = TRUE,
  use_log_precip_on_wet = TRUE,
  method = "pearson",
  min_pairs = 50
)
```

## Arguments

- data:

  Data frame with daily values and \`date\`/\`id\` columns.

- variables:

  Character vector of variable names (must include \`"precip"\`).

- mc_thresholds:

  Optional wet-day thresholds (required for \`wet_def =
  "monthly_quantile"\`).

- wet_def:

  Wet-day definition: \`"gt0"\` or \`"monthly_quantile"\`.

- use_anom:

  Logical; subtract monthly means before correlation.

- use_log_precip_on_wet:

  Logical; apply \`log1p\` to precip on wet days.

- method:

  Correlation method passed to \`stats::cor\`.

- min_pairs:

  Minimum number of paired values required to report a correlation.

## Value

Data frame of conditional correlations by grid, regime, and variable.
