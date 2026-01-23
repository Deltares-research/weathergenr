# Compute recommended maximum MODWT levels based on length and filter width

Compute recommended maximum MODWT levels based on length and filter
width

## Usage

``` r
.recommended_max_levels(n, filter_length, max_period_frac = 1)
```

## Arguments

- n:

  Integer scalar. Series length.

- filter_length:

  Integer scalar. Wavelet filter length.

- max_period_frac:

  Numeric scalar in (0, 1\]. Optional cap based on record-length
  fraction.

## Value

Integer scalar. Recommended maximum levels (J), at least 1.
