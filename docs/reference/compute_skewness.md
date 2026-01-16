# Compute sample skewness

Computes the (non-bias-corrected) third standardized central moment
(skewness) for a numeric vector, excluding `NA` values.

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
