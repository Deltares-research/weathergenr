# Vectorized variance matching for simulation matrices

If the relative standard deviation mismatch exceeds tol, rescales each
column to match target_sd while keeping the column mean at target_mean.

## Usage

``` r
.variance_match_matrix(sim, target_sd, tol, target_mean)
```

## Arguments

- sim:

  Numeric matrix.

- target_sd:

  Numeric scalar.

- tol:

  Numeric scalar in \[0, 1\].

- target_mean:

  Numeric scalar.

## Value

Numeric matrix with the same dimension as sim.
