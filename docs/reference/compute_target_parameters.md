# Compute Target Gamma Parameters from Baseline and Change Factors

Computes year-specific target gamma distribution parameters for each
month by scaling baseline monthly mean and variance using change
factors, then converting target moments to gamma \`shape\` and
\`scale\`.

## Usage

``` r
compute_target_parameters(
  base_gamma,
  mean_factor,
  var_factor,
  months_ok,
  n_years
)
```

## Arguments

- base_gamma:

  Data frame produced by \`fit_monthly_distributions()\`, with columns
  \`month\`, \`mean\`, and \`var\` (and typically \`shape\`, \`scale\`).

- mean_factor:

  Numeric matrix of multiplicative mean factors with dimensions
  \`n_years x 12\` (year x month).

- var_factor:

  Numeric matrix of multiplicative variance factors with dimensions
  \`n_years x 12\` (year x month).

- months_ok:

  Integer vector of months expected to be represented (kept for
  interface consistency; the function uses \`base_gamma\$month\`).

- n_years:

  Integer. Number of years (must match \`nrow(mean_factor)\` and
  \`nrow(var_factor)\`).

## Value

A list with elements:

- \`months\` Integer vector of months (from \`base_gamma\$month\`).

- \`shape\` Numeric matrix \`\[n_months x n_years\]\` of target shapes.

- \`scale\` Numeric matrix \`\[n_months x n_years\]\` of target scales.

- \`mean\` Numeric matrix \`\[n_months x n_years\]\` of target means.

- \`var\` Numeric matrix \`\[n_months x n_years\]\` of target variances.
