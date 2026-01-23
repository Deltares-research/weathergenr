# Block bootstrap for many realizations

Generates n_sim bootstrap realizations using contiguous blocks sampled
with replacement. The implementation is vectorized over realizations
within each block.

## Usage

``` r
.block_bootstrap_matrix(x, n, n_sim, block_len)
```

## Arguments

- x:

  Numeric vector.

- n:

  Integer. Output length.

- n_sim:

  Integer. Number of realizations.

- block_len:

  Integer. Block length.

## Value

Numeric matrix n by n_sim.
