# Log initial setup information

Internal helper. Displays general information at the start of filtering.

## Usage

``` r
log_filtering_start(
  n_obs,
  n_sim,
  n_realizations,
  sample_target,
  relax_priority
)
```

## Arguments

- n_obs:

  Integer scalar. Number of observations in observed series.

- n_sim:

  Integer scalar. Number of years in simulated series used after
  windowing.

- n_realizations:

  Integer scalar. Number of candidate realizations.

- sample_target:

  Integer scalar. Target number to select.

- relax_priority:

  Character vector. Relaxation priority vector.

## Value

Invisibly returns NULL.
