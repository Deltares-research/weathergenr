# Fit Monthly Gamma Distributions for Wet-Day Precipitation

Fits a gamma distribution to wet-day precipitation amounts for each
month in \`months_ok\`. Uses \`fitdistrplus::fitdist()\` and returns
month-specific gamma parameters and implied moments. Months that fail to
fit are dropped.

## Usage

``` r
fit_monthly_distributions(
  precip_wet,
  month_wet,
  months_ok,
  fit_method,
  verbose
)
```

## Arguments

- precip_wet:

  Numeric vector of wet-day precipitation values (strictly \> 0
  recommended). Length must match \`month_wet\`.

- month_wet:

  Integer vector of months (1-12) corresponding to each element of
  \`precip_wet\`.

- months_ok:

  Integer vector of months to attempt fitting (subset of 1-12).

- fit_method:

  Character scalar. Estimation method forwarded to
  \`fitdistrplus::fitdist()\` (e.g., \`"mle"\`, \`"mme"\`).

- verbose:

  Logical. If \`TRUE\`, emits a warning when a monthly fit fails.

## Value

A \`data.frame\` with columns:

- \`month\` Month number.

- \`shape\` Gamma shape parameter.

- \`scale\` Gamma scale parameter.

- \`mean\` Implied mean (\`shape \* scale\`).

- \`var\` Implied variance (\`shape \* scale^2\`).

Months with failed fits are removed. The returned data frame has an
attribute \`n_failed_fits\` with the number of failed months.
