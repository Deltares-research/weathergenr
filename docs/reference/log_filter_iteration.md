# Log one filtering iteration summary

Writes a compact per iteration summary showing pass counts, pass rates,
and current criteria for each active filter family, followed by the
current pool size relative to the target.

## Usage

``` r
log_filter_iteration(
  iter,
  passes,
  pool,
  n_total,
  target,
  bounds,
  tail_metrics,
  wavelet_active,
  spectral_diag,
  note = NULL
)
```

## Arguments

- iter:

  Integer scalar. Iteration number. Use 0 for the initial evaluation.

- passes:

  Named list. Logical vectors indicating pass or fail for each filter
  family.

- pool:

  Integer vector. Indices of realizations currently in the pool.

- n_total:

  Integer scalar. Total number of realizations evaluated.

- target:

  Integer scalar. Target pool size.

- bounds:

  Environment. Current bounds values.

- tail_metrics:

  List. Tail metrics produced by
  [`compute_tailmass_metrics()`](https://deltares-research.github.io/weathergenr/reference/compute_tailmass_metrics.md).

- wavelet_active:

  Logical scalar. TRUE if wavelet filtering is active.

- spectral_diag:

  List. Spectral diagnostics returned by
  [`compute_spectral_metrics()`](https://deltares-research.github.io/weathergenr/reference/compute_spectral_metrics.md).

- note:

  Character scalar or NULL. Optional message describing the action taken
  in this iteration.

## Value

Invisibly returns NULL.
