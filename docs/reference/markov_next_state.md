# Determine Next State in a First-Order Markov Chain

Given the previous state and a random number \`u_rand\`, this function
returns the next state in a 3-state first-order Markov chain. The
transition probabilities are state- and index-specific.

## Usage

``` r
markov_next_state(state_prev, u_rand, idx, p00, p01, p10, p11, p20, p21)
```

## Arguments

- state_prev:

  Integer. The current (previous) state (0, 1, or 2).

- u_rand:

  Numeric. A random number in \[0, 1\] used to sample the next state.

- idx:

  Integer. Index for selecting transition probabilities (e.g., time step
  or spatial location).

- p00:

  Numeric vector. Probability of transition from state 0 to 0.

- p01:

  Numeric vector. Probability of transition from state 0 to 1.

- p10:

  Numeric vector. Probability of transition from state 1 to 0.

- p11:

  Numeric vector. Probability of transition from state 1 to 1.

- p20:

  Numeric vector. Probability of transition from state 2 to 0.

- p21:

  Numeric vector. Probability of transition from state 2 to 1.

## Value

Integer. The next state (0, 1, or 2).

## Details

States are typically defined as: - 0: Dry - 1: Wet - 2: Extreme

The transition is determined using cumulative probabilities for the
transitions to state 0 and 1. Any remaining probability is assigned to
state 2.

## Examples

``` r
set.seed(123)
u_rand <- runif(1)
markov_next_state(
  state_prev = 1,
  u_rand = u_rand,
  idx = 5,
  p00 = rep(0.7, 10),
  p01 = rep(0.2, 10),
  p10 = rep(0.3, 10),
  p11 = rep(0.4, 10),
  p20 = rep(0.1, 10),
  p21 = rep(0.3, 10)
)
#> [1] 0
```
