# Safely Get or Sample an Index from a Vector

Returns a valid index from \`candidate_prcp\` based on the provided
\`idx\`. If \`idx\` is missing, out of bounds, or invalid, a random
index is sampled instead. If \`candidate_prcp\` is empty,
\`NA_integer\_\` is returned.

## Usage

``` r
get_result_index(idx, candidate_prcp)
```

## Arguments

- idx:

  Integer. Proposed index for retrieving an element from
  \`candidate_prcp\`.

- candidate_prcp:

  Numeric vector. A set of candidate precipitation values (e.g., for the
  next day).

## Value

Integer. A valid index from \`candidate_prcp\`, or \`NA_integer\_\` if
\`candidate_prcp\` is empty.

## Examples

``` r
candidate_prcp <- c(0, 5, 10, 20)
get_result_index(2, candidate_prcp)
#> [1] 2

set.seed(1)
get_result_index(10, candidate_prcp)
#> [1] 1

get_result_index(NA, candidate_prcp)
#> [1] 4
get_result_index(1, numeric(0))
#> [1] NA
```
