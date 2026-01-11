# Relax bounds for one filter

Applies one relaxation step to a single filter. Updates bounds in place
and returns status.

## Usage

``` r
relax_bounds_one_filter(
  filter_name,
  bounds_env,
  wavelet_active_env,
  recompute_tailmass_fn
)
```

## Arguments

- filter_name:

  Character. Name of filter to relax

- bounds_env:

  Environment containing bounds

- wavelet_active_env:

  Environment containing wavelet_active flag

- recompute_tailmass_fn:

  Function to recompute tail mass when thresholds change

## Value

List with changed (logical) and msg (character)
