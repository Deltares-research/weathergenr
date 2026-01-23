# Log filtering setup information

Writes a header block describing the filtering configuration. This
includes series lengths, the number of candidate realizations, the
selection target, and the relaxation priority.

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

  Integer scalar. Length of the observed series in years.

- n_sim:

  Integer scalar. Length of the simulated series in years.

- n_realizations:

  Integer scalar. Number of candidate realizations.

- sample_target:

  Integer scalar. Number of realizations requested.

- relax_priority:

  Character vector. Relaxation priority ordering.

## Value

Invisibly returns NULL.
