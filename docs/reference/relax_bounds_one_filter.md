# Relax bounds for one filter family

Applies one relaxation step for a single filter family and updates the
bounds environment in place. This function is called by
[`filter_warm_pool()`](https://deltares-research.github.io/weathergenr/reference/filter_warm_pool.md).

Relaxation behavior by filter family:

- mean:

  Increases the mean tolerance up to `relax_mean_max`.

- sd:

  Increases the sd tolerance up to `relax_sd_max`.

- tail_low:

  First increases tail tolerance, then increases `tail_low_p` in steps,
  recomputing tail metrics after changes.

- tail_high:

  First increases tail tolerance, then decreases `tail_high_p` in steps,
  recomputing tail metrics after changes.

- wavelet:

  Relaxes spectral correlation threshold, then relaxes peak match
  fraction, then disables peak matching, then disables wavelet
  filtering.

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

  Character scalar. Filter family to relax. Must be one of `"mean"`,
  `"sd"`, `"tail_low"`, `"tail_high"`, `"wavelet"`.

- bounds_env:

  Environment. Environment containing the current bounds using snake
  case keys. Updated in place.

- wavelet_active_env:

  Environment. Environment containing a logical `wavelet_active` flag.
  May be updated in place.

- recompute_tailmass_fn:

  Function. Callback used to recompute tail mass metrics when tail
  quantile thresholds change.

## Value

Named list with: `changed` logical scalar indicating whether a bound was
changed, `msg` character message describing the applied change.
