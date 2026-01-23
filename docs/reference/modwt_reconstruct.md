# Inverse MODWT

Reconstructs a time series from MODWT coefficients.

## Usage

``` r
modwt_reconstruct(wt)
```

## Arguments

- wt:

  Either a `"modwt_result"` object returned by
  [`modwt_decompose`](https://deltares-research.github.io/weathergenr/reference/modwt_decompose.md),
  or a raw
  [`waveslim::modwt()`](https://rdrr.io/pkg/waveslim/man/modwt.html)
  output list containing W\* and VJ (or V\<level\>) entries.

## Value

Numeric vector. Reconstructed time series.
