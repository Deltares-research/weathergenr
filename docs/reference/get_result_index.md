# Safely Get or Sample an Index from a Vector

Returns a valid index from \`candidate_precip\` based on the provided
\`idx\`. If \`idx\` is missing, out of bounds, or invalid, a random
index is sampled instead. If \`candidate_precip\` is empty,
\`NA_integer\_\` is returned.

## Usage

``` r
get_result_index(idx, candidate_precip)
```

## Arguments

- idx:

  Integer-valued scalar. Proposed index for retrieving an element from
  \`candidate_precip\`. Non-integer or invalid values trigger
  resampling.

- candidate_precip:

  Numeric vector. A set of candidate precipitation values (e.g., for the
  next day).

## Value

Integer. A valid index from \`candidate_precip\`, or \`NA_integer\_\` if
\`candidate_precip\` is empty.

## Examples

``` r
candidate_precip <- c(0, 5, 10, 20)
get_result_index(2, candidate_precip)
#> [1] 2

set.seed(1)
get_result_index(10, candidate_precip)
#> [1] 1

get_result_index(NA, candidate_precip)
#> [1] 4
get_result_index(1, numeric(0))
#> [1] NA
```
