# Log filter iteration details

Internal helper. Prints iteration diagnostics in a compact table-like
format.

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
  wavelet_pars,
  note = NULL
)
```

## Arguments

- iter:

  Integer scalar. Iteration number.

- passes:

  Named list of logical vectors. Per-filter pass vectors.

- pool:

  Integer vector. Pool indices passing all active filters.

- n_total:

  Integer scalar. Total number of realizations.

- target:

  Integer scalar. Target pool size.

- bounds:

  Environment or list of bounds.

- tail_metrics:

  List. Tail metrics used for criteria display.

- wavelet_active:

  Logical scalar.

- wavelet_pars:

  List. Wavelet parameter list.

- note:

  Optional character scalar.

## Value

Invisibly returns NULL.
