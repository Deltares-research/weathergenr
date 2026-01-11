# Quantile Mapping for Climate Change Perturbation of Precipitation

Applies quantile mapping to adjust a daily precipitation time series to
reflect changes in monthly mean and variance, consistent with climate
change scenarios. The method fits a gamma distribution to observed
nonzero precipitation in each calendar month, then modifies the
distribution's mean and variance according to supplied monthly/yearly
change factors, and maps the original values to their new quantiles in
the perturbed distribution.

This function is designed for climate stress-testing workflows to impose
user-specified changes in mean and variance on daily precipitation data
while preserving realistic distributional characteristics.

## Usage

``` r
quantile_mapping(
  value = NULL,
  mean.change = NULL,
  var.change = NULL,
  mon.ts = NULL,
  year.ts = NULL,
  fit.method = "mme",
  min.events = 10,
  validate.output = TRUE,
  verbose = FALSE,
  compute.diagnostics = FALSE
)
```

## Arguments

- value:

  Numeric vector. Original daily precipitation values to perturb
  (typically for one grid cell). Must be non-negative.

- mean.change:

  Numeric matrix of mean change factors (multiplicative), dimension:
  \`n_years\` x 12 (year, month). Each entry scales that year-month's
  mean.

- var.change:

  Numeric matrix of variance change factors (multiplicative), dimension:
  \`n_years\` x 12 (year, month). Each entry scales that year-month's
  variance.

- mon.ts:

  Integer vector (same length as \`value\`). Calendar month for each day
  (1-12).

- year.ts:

  Integer vector (same length as \`value\`). Simulation year index for
  each day (1 = first year, etc).

- fit.method:

  Character. Method for fitting the base gamma distribution; passed to
  \[fitdistrplus::fitdist()\]. Default is \`"mme"\` (method of moments).

- min.events:

  Integer. Minimum number of non-zero precipitation events required per
  month for distribution fitting. Default is 10.

- validate.output:

  Logical. If TRUE, checks output for NaN/Inf and replaces with original
  values. Default is TRUE.

- verbose:

  Logical. If TRUE, prints diagnostic messages about months that cannot
  be perturbed. Default is FALSE.

## Value

Numeric vector, same length as \`value\`. Precipitation time series
perturbed according to quantile mapping procedure. Includes attributes:

- \`perturbed_months\`: Integer vector of months that were successfully
  perturbed

- \`skipped_months\`: Integer vector of months skipped due to
  insufficient data

- \`n_failed_fits\`: Number of distribution fits that failed

## Details

\## Distribution Fitting The function fits a gamma distribution to
non-zero precipitation in each calendar month. Months with fewer than
\`min.events\` non-zero days are skipped. The gamma distribution is
parameterized with shape and scale parameters derived from the method of
moments or maximum likelihood estimation.

\## Quantile Mapping Procedure 1. Fit base gamma distribution to
observed non-zero precipitation by month 2. Compute target distribution
parameters by scaling base mean and variance 3. Map each original value
to its quantile in the base distribution 4. Transform to the same
quantile in the target distribution

\## Limitations - Zero precipitation values remain zero (dry day
frequency unchanged) - Requires sufficient non-zero events per month for
reliable fitting - Extrapolation beyond observed range may produce
unrealistic extreme values - Change factors must be positive; negative
or zero values will cause errors

## See also

[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example: 2 years of daily data, 5% mean and 10% variance increase
set.seed(123)
n_days <- 730
year_idx <- rep(1:2, each = 365)
month_idx <- rep(rep(1:12, times = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)), 2)[1:n_days]
daily_precip <- rgamma(n_days, shape = 1, scale = 5)

mean.change <- matrix(1.05, nrow = 2, ncol = 12)
var.change <- matrix(1.10, nrow = 2, ncol = 12)

perturbed_precip <- quantile_mapping(
  value = daily_precip,
  mean.change = mean.change,
  var.change = var.change,
  mon.ts = month_idx,
  year.ts = year_idx
)

# Check which months were perturbed
attr(perturbed_precip, "perturbed_months")
attr(perturbed_precip, "skipped_months")
} # }
```
