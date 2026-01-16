# Compute sample excess kurtosis

Computes the (non-bias-corrected) fourth standardized central moment
minus 3 (excess kurtosis) for a numeric vector, excluding `NA` values.

## Usage

``` r
compute_kurtosis(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric scalar. Excess kurtosis estimate, or `NA`.

## Details

Returns `NA` if fewer than 4 non-missing values are available or if the
standard deviation is zero.
