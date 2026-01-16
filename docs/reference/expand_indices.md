# Expand Index Positions by Fixed Offsets

Expands a vector of base index positions by adding a set of integer
offsets, returning all valid index positions within a specified upper
bound. This is typically used to generate moving windows or neighborhood
indices around reference positions (e.g., expanding event days to
adjacent days).

## Usage

``` r
expand_indices(base_idx, offset_vec, n_max)
```

## Arguments

- base_idx:

  Integer vector of base index positions to be expanded.

- offset_vec:

  Integer vector of offsets to apply to each element of `base_idx`.
  Offsets may be positive, negative, or zero.

- n_max:

  Integer. Maximum allowable index length; expanded indices are
  constrained to the range `1:(n_max - 1)`.

## Value

Integer vector of expanded index positions that satisfy the validity
constraints.

## Details

For each element in `base_idx`, all values in `offset_vec` are added.
Resulting indices that are less than or equal to zero, or for which
`index + 1 > n_max`, are discarded. The function does not enforce
uniqueness or sorting of the output.

## Examples

``` r
base_idx <- c(5, 10)
offset_vec <- -1:1
weathergenr:::expand_indices(base_idx, offset_vec, n_max = 15)
#> [1]  4  5  6  9 10 11
```
