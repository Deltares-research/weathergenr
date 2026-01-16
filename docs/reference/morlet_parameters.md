# Compute Morlet Wavelet Parameters

Computes key parameters for Morlet wavelets: Fourier factor
(scale-\>period), cone-of-influence (COI) e-folding time, and minimum
degrees of freedom.

## Usage

``` r
morlet_parameters(k0 = 6)
```

## Arguments

- k0:

  Numeric scalar. Morlet nondimensional frequency (default: 6).

## Value

Named numeric vector: `fourier_factor`, `coi`, `dofmin`.
