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
