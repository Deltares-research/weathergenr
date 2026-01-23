# Fast ARMA simulation for many realizations

Simulates a stationary ARMA process using direct recursion with Gaussian
innovations. This is intended for small p and q (0 to 2) and large
n_sim.

## Usage

``` r
.simulate_arma_matrix(ar, ma, sd, n, n_sim, burnin = 100L)
```

## Arguments

- ar:

  Numeric vector. AR coefficients.

- ma:

  Numeric vector. MA coefficients.

- sd:

  Numeric scalar. Innovation standard deviation.

- n:

  Integer. Output length.

- n_sim:

  Integer. Number of realizations.

- burnin:

  Integer. Burn in length.

## Value

Numeric matrix n by n_sim.
