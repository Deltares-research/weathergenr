# Log filter iteration details with table format

Prints iteration diagnostics in table format for all iterations.

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

  Iteration number

- passes:

  List of pass vectors

- pool:

  Vector of pool indices

- n_total:

  Total number of realizations

- target:

  Target pool size

- bounds:

  Bounds environment or list

- tail_metrics:

  Tail metrics list

- wavelet_active:

  Logical

- wavelet_pars:

  Wavelet parameters list

- note:

  Optional note string
