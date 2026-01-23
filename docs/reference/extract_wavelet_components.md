# Reconstruct Wavelet Components from Selected Scales (CWT)

Reconstructs time-domain components from selected wavelet scales using
the Torrence & Compo (1998) inverse transform approximation. Optionally
adds a residual component representing non-selected scales to ensure
additive closure.

Note: CWT-reconstructed components are additive but NOT orthogonal due
to scale overlap. For MODWT additive decomposition, use
[`analyze_wavelet_additive`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_additive.md).

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

  Numeric scalar. psi_0(0) constant (default: pi^(-1/4)).

- include_residual:

  Logical. If `TRUE`, adds a residual (non-selected) component.

## Value

Numeric matrix with one column per selected component and, if requested,
a final `"Noise"` residual column.

## References

Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61-78.

## See also

[`analyze_wavelet_spectrum`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md),
[`analyze_wavelet_additive`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_additive.md)
