# Regrid Global Wavelet Spectrum to a Target Period Grid

Linearly interpolates the global wavelet spectrum from a wavelet
analysis output onto a target period grid and fills edge values using
nearest-neighbor filling.

## Usage

``` r
gws_regrid(wavelet, target_period, use_unmasked = FALSE)
```

## Arguments

- wavelet:

  List output from
  [`analyze_wavelet_spectrum`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md).

- target_period:

  Numeric vector. Target periods for interpolation.

- use_unmasked:

  Logical. If `TRUE` and available, uses `wavelet$gws_unmasked`;
  otherwise uses `wavelet$gws`.

## Value

Numeric vector of interpolated GWS values on `target_period`.
