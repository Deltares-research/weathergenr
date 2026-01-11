# Safely Get or Sample an Index from a Vector

Returns a valid index from \`precip.tomorrow\` based on the provided
\`result\`. If \`result\` is missing, out of bounds, or invalid, a
random index is sampled instead. If \`precip.tomorrow\` is empty,
\`NA_integer\_\` is returned.

## Usage

``` r
get_result_index(result, precip.tomorrow)
```

## Arguments

- result:

  Integer. Proposed index for retrieving an element from
  \`precip.tomorrow\`.

- precip.tomorrow:

  Numeric vector. A set of candidate precipitation values (e.g., for the
  next day).

## Value

Integer. A valid index from \`precip.tomorrow\`, or \`NA_integer\_\` if
\`precip.tomorrow\` is empty.

## Details

This function is useful in stochastic weather generation workflows where
fallback behavior is needed when a computed index is not valid.

## Examples

``` r
precip.tomorrow <- c(0, 5, 10, 20)

# Valid index
get_result_index(2, precip.tomorrow) # returns 2
#> [1] 2

# Invalid index: falls back to sampling
set.seed(1)
get_result_index(10, precip.tomorrow) # randomly returns a valid index
#> [1] 1

# NA index: falls back to sampling
get_result_index(NA, precip.tomorrow)
#> [1] 4

# Empty vector: returns NA
get_result_index(1, numeric(0))
#> [1] NA
```
