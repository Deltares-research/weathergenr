# Apply Gamma-to-Gamma Quantile Mapping for Precipitation

Internal helper implementing gamma-based quantile mapping for
precipitation: values are mapped to probability space under a baseline
monthly gamma, then mapped back to precipitation under a target
monthly/year gamma. Optionally applies a tail-probability exaggeration
and an enforcement step to match the target mean at the month-year
subset level.

## Usage

``` r
apply_quantile_mapping(
  precip,
  mon,
  year,
  base_gamma,
  target_gamma,
  exaggerate_extremes = FALSE,
  extreme_prob_threshold = 0.95,
  extreme_k = 1.2,
  enforce_target_mean = TRUE
)
```

## Arguments

- precip:

  Numeric vector of precipitation values (subset, typically wet days).

- mon:

  Integer vector of months (1-12), same length as \`precip\`.

- year:

  Integer vector of year indices used to select columns in
  \`target_gamma\$shape/scale/mean\`. This is an index (1..n_years), not
  necessarily the calendar year value.

- base_gamma:

  Data frame with \`month\`, \`shape\`, and \`scale\` columns.

- target_gamma:

  List produced by \`compute_target_parameters()\`, containing
  \`shape\`, \`scale\`, and \`mean\` matrices indexed as \`\[month_row,
  year_index\]\`.

- exaggerate_extremes:

  Logical. If \`TRUE\`, modifies probabilities above
  \`extreme_prob_threshold\` as \`u' = 1 - (1-u)^k\`.

- extreme_prob_threshold:

  Numeric in (0,1). Tail threshold \`u0\`.

- extreme_k:

  Positive numeric. Tail exponent \`k\`.

- enforce_target_mean:

  Logical. If \`TRUE\`, rescales mapped values within each (month, year)
  subset to match \`target_gamma\$mean\`.

## Value

Numeric vector of mapped precipitation values, same length as
\`precip\`.
