# Filter correlations by type (cross-grid or inter-variable)

Filter correlations by type (cross-grid or inter-variable)

## Usage

``` r
.filter_correlation_pairs(
  cor_data,
  same_var,
  same_id,
  allowed_pairs,
  pair_type
)
```

## Arguments

- cor_data:

  Data frame of correlations.

- same_var:

  Logical; require same-variable pairs.

- same_id:

  Logical; require same-grid pairs.

- allowed_pairs:

  Character vector of allowed pair keys.

- pair_type:

  \`"id"\` for grid pairs or \`"variable"\` for variable pairs.

## Value

Filtered data frame of correlations.
