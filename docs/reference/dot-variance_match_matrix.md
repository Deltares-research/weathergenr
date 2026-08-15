# Vectorized variance matching for simulation matrices

Where a column's sample standard deviation differs from \`target_sd\` by
more than \`tol\`, rescales that column about \*\*its own mean\*\* so
the deviation is corrected and the mean is left alone.

## Usage

``` r
.variance_match_matrix(sim, target_sd, tol)
```

## Arguments

- sim:

  Numeric matrix, one realization per column.

- target_sd:

  Numeric scalar. Sample standard deviation to match.

- tol:

  Numeric scalar in \[0, 1\]. Relative tolerance below which a column is
  left untouched.

## Value

Numeric matrix with the same dimension as sim.

## Details

This function previously re-centred each corrected column on a scalar
\`target_mean\`, which replaced the column's mean rather than preserving
it. Because roughly half an ensemble typically exceeds \`tol\`, that put
a point mass on the observed mean: those realizations carried no spread
in their own mean at all, and \`filter_warm_pool()\`'s mean criterion
could not reject any of them.

That spread is not noise to be removed. The sampling distribution of a
20-40 year mean under a persistent process is precisely the
low-frequency variability WARM exists to reproduce, so collapsing it
works against the method. Rescaling about the column's own mean is also
what "variance matching" says on the tin, and it makes the corrected and
uncorrected columns agree about how the mean is carried – previously
they did not.
