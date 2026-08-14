# Compute sample excess kurtosis

Computes sample excess kurtosis as \\\sum (x-\bar{x})^4 / (n s^4) - 3\\
for a numeric vector, excluding `NA` values.

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
standard deviation is zero. As in
[`compute_skewness()`](https://deltares-research.github.io/weathergenr/reference/compute_skewness.md),
the standardization uses
[`stats::sd()`](https://rdrr.io/r/stats/sd.html) (denominator \\n-1\\).
