# Morlet Wavelet in Fourier Domain

Computes the Morlet wavelet "daughter" in Fourier space for a given
scale. Normalization follows Torrence & Compo (1998).

## Usage

``` r
morlet_wavelet(k, scale, k0 = 6)
```

## Arguments

- k:

  Numeric vector. Angular frequencies in FFT ordering.

- scale:

  Numeric scalar. Wavelet scale (inverse frequency), must be positive.

- k0:

  Numeric scalar. Morlet nondimensional frequency (default: 6).

## Value

Numeric vector (complex). Morlet wavelet in Fourier space.

## References

Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61-78.

## See also

[`morlet_parameters`](https://deltares-research.github.io/weathergenr/reference/morlet_parameters.md),
[`analyze_wavelet_spectrum`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md)
