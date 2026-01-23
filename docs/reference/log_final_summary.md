# Log final filtering summary

Writes a footer block summarizing the final pool size, number sampled,
and the selection mode that was used.

## Usage

``` r
log_final_summary(pool_size, n_total, n_sampled, relaxation_level)
```

## Arguments

- pool_size:

  Integer scalar. Final candidate pool size.

- n_total:

  Integer scalar. Total number of candidate realizations.

- n_sampled:

  Integer scalar. Number of realizations sampled.

- relaxation_level:

  Character scalar. Label describing the selection mode.

## Value

Invisibly returns NULL.
