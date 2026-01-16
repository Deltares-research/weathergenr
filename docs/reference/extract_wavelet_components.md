# Reconstruct Wavelet Components from Selected Scales

Reconstructs time-domain components from selected wavelet scales using
the Torrence & Compo (1998) inverse transform approximation. Optionally
adds a residual component representing non-selected scales to ensure
additive closure.

## Usage

``` r
extract_wavelet_components(
  wave,
  signif_periods,
  scale,
  dj = 0.25,
  dt = 1,
  series_sd,
  series_mean = 0,
  Cdelta = 0.776,
  w0_0 = pi^(-1/4),
  include_residual = TRUE
)
```

## Arguments

- wave:

  Complex matrix (scales x time). Wavelet coefficients.

- signif_periods:

  Integer vector. Scale indices to reconstruct.

- scale:

  Numeric vector. Scale values aligned with `nrow(wave)`.

- dj:

  Numeric scalar. Scale resolution used in the transform.

- dt:

  Numeric scalar. Time step of the series.

- series_sd:

  Numeric scalar. Standard deviation of the original series
  (pre-standardization).

- series_mean:

  Numeric scalar. Mean of the original series; added back to residual.

- Cdelta:

  Numeric scalar. Morlet reconstruction constant (default: 0.776).

- w0_0:

  Numeric scalar. \\\psi_0(0)\\ constant (default: \\\pi^{-1/4}\\).

- include_residual:

  Logical. If `TRUE`, adds a residual (non-selected) component.

## Value

Numeric matrix with one column per selected component and, if requested,
a final `"Noise"` residual column.
