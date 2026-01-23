# Compute max allowed MODWT levels given a max-period fraction

Ensures the smooth threshold period 2^(J+1) does not exceed floor(n \*
max_period_frac). This prevents representing variability at periods
longer than a chosen fraction of the record length.

## Usage

``` r
.max_levels_by_period_frac(n, max_period_frac)
```

## Arguments

- n:

  Integer scalar. Series length.

- max_period_frac:

  Numeric scalar in (0, 1\]. Maximum allowed period as a fraction of n.

## Value

Integer scalar. Maximum allowed number of levels (J), at least 1.
