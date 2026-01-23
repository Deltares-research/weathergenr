# Choose n_levels given requested and recommended maximum

Choose n_levels given requested and recommended maximum

## Usage

``` r
.select_n_levels(n_levels, max_levels, default_backoff = 1L)
```

## Arguments

- n_levels:

  Integer scalar or NULL.

- max_levels:

  Integer scalar.

- default_backoff:

  Integer scalar. If NULL n_levels, use max(1, max_levels -
  default_backoff).

## Value

Integer scalar. Selected n_levels in \[1, max_levels\].
