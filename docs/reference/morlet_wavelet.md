# Morlet Wavelet in Fourier Domain

Computes the Morlet wavelet function in Fourier space for a given scale.
This is the "daughter" wavelet used in continuous wavelet transform
analysis. The function applies normalization following Torrence & Compo
(1998).

## Usage

``` r
morlet_wavelet(k, s, k0 = 6)
```

## Arguments

- k:

  Numeric vector. Wave number vector (angular frequency).

- s:

  Numeric scalar. Scale parameter (inverse of frequency).

- k0:

  Numeric scalar. Omega0 parameter for the Morlet wavelet (default = 6).
  This value balances time and frequency localization.

## Value

Numeric vector (complex). The Morlet wavelet in Fourier space, with
length equal to length(k).

## Details

The Morlet wavelet is defined in Fourier space as: psi(s k) = pi^(-1/4)
\* sqrt(s k_2 n) \* exp(-(s k - k_0)^2 / 2) where the exponential term
is applied only to positive frequencies.

The normalization factor ensures energy conservation and proper inverse
transform reconstruction.

## References

Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet
Analysis. Bulletin of the American Meteorological Society, 79(1), 61-78.
