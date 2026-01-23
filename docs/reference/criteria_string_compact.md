# Create compact criteria text for logging

Builds a compact one line string describing the current criteria for a
given filter family. Used by the logging helpers in this script.

## Usage

``` r
criteria_string_compact(
  filter_name,
  bounds,
  tail_metrics,
  wavelet_active,
  spectral_diag
)
```

## Arguments

- filter_name:

  Character scalar. Filter family name.

- bounds:

  List or environment. Bounds values used to build the criteria string.

- tail_metrics:

  List. Tail metrics produced by
  [`compute_tailmass_metrics()`](https://deltares-research.github.io/weathergenr/reference/compute_tailmass_metrics.md).

- wavelet_active:

  Logical scalar. TRUE if wavelet filtering is active.

- spectral_diag:

  List. Spectral diagnostics returned by
  [`compute_spectral_metrics()`](https://deltares-research.github.io/weathergenr/reference/compute_spectral_metrics.md).

## Value

Character scalar. A compact criteria string.
