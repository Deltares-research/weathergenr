# Compute sample skewness

Computes sample skewness as \\\sum (x-\bar{x})^3 / (n s^3)\\ for a
numeric vector, excluding `NA` values.

## Usage

``` r
compute_skewness(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric scalar. Skewness estimate, or `NA`.

## Details

Returns `NA` if fewer than 3 non-missing values are available or if the
standard deviation is zero (which yields undefined standardization).
This estimator uses [`stats::sd()`](https://rdrr.io/r/stats/sd.html)
(denominator \\n-1\\) in the standardization term and follows the same
convention used by `e1071::skewness(type = 2)`.
