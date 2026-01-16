# Normalize a Probability Vector with Robust Fallback Handling

Normalizes a numeric vector to sum to one, treating non-finite and
negative values as zero. If the total mass is zero after cleaning, a
fallback distribution is returned instead.

## Usage

``` r
normalize_probs(prob, fallback_prob = NULL)
```

## Arguments

- prob:

  Numeric vector of (unnormalized) probabilities or weights. Non-finite
  values (`NA`, `NaN`, `Inf`) and negative values are set to zero prior
  to normalization.

- fallback_prob:

  Optional numeric vector to return when `prob` has zero total mass
  after cleaning. If `NULL` (default), a uniform distribution of length
  `length(prob)` is returned.

## Value

Numeric vector of probabilities summing to one. If a fallback is used,
the returned vector is either `fallback_prob` or a uniform distribution.

## Details

This helper is intended for internal use in stochastic or probabilistic
workflows where degenerate or numerically unstable probability vectors
may arise. The function does not check that `fallback_prob` sums to one;
this is assumed to be ensured by the caller.

## Examples

``` r
if (FALSE) { # \dontrun{
prob <- c(0.2, NA, -0.1, 0.5)
normalize_probs(prob)

normalize_probs(c(0, 0, 0))

normalize_probs(c(0, 0, 0), fallback_prob = c(0.7, 0.2, 0.1))
} # }
```
