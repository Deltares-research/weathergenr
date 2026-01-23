# Compute Morlet Wavelet Parameters

Computes key parameters for Morlet wavelets: Fourier factor (scale to
period), cone-of-influence (COI) e-folding time, and minimum degrees of
freedom.

## Usage

``` r
morlet_parameters(k0 = 6)
```

## Arguments

- k0:

  Numeric scalar. Morlet nondimensional frequency (default: 6).

## Value

Named numeric vector: `fourier_factor`, `coi`, `dofmin`.

## References

Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61-78.

## See also

[`morlet_wavelet`](https://deltares-research.github.io/weathergenr/reference/morlet_wavelet.md),
[`analyze_wavelet_spectrum`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md)
