# Compute Validation Diagnostics for Quantile Mapping

Calculates comprehensive validation metrics comparing original and
adjusted precipitation time series. Evaluates moment preservation,
distributional characteristics, extreme values, and temporal patterns.

## Usage

``` r
validate_quantile_mapping(
  value.original,
  value.adjusted,
  mon.ts = NULL,
  year.ts = NULL,
  mean.change = NULL,
  var.change = NULL,
  wet.threshold = 0.1,
  quantiles = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
)
```

## Arguments

- value.original:

  Numeric vector. Original precipitation values.

- value.adjusted:

  Numeric vector. Adjusted precipitation values (same length).

- mon.ts:

  Integer vector. Month indices (1-12) for each value.

- year.ts:

  Integer vector. Year indices for each value.

- mean.change:

  Numeric matrix. Target mean change factors (n_years x 12).

- var.change:

  Numeric matrix. Target variance change factors (n_years x 12).

- wet.threshold:

  Numeric. Threshold for defining wet days (default = 0.1 mm).

- quantiles:

  Numeric vector. Quantiles to evaluate (default = standard set).

## Value

A list of class "qmap_diagnostics" containing:

- `moments`: Data frame of moment statistics

- `quantiles`: Data frame of quantile comparisons

- `extremes`: Data frame of extreme value metrics

- `temporal`: Data frame of temporal pattern metrics

- `monthly`: List of monthly-specific diagnostics

- `spells`: Data frame of spell length statistics

- `summary`: Overall assessment metrics
