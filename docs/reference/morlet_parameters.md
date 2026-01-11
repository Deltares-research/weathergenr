# Compute Morlet Wavelet Parameters

Calculates key parameters for the Morlet wavelet transform, including
the Fourier factor (scale-to-period conversion), cone of influence
(COI), and minimum degrees of freedom.

## Usage

``` r
morlet_parameters(k0 = 6)
```

## Arguments

- k0:

  Numeric scalar. Omega0 parameter for the Morlet wavelet (default = 6).
  Standard choice is 6, which gives ~6 oscillations per wavelet.

## Value

Named numeric vector with three elements: fourier_factor: Conversion
factor from wavelet scale to Fourier period coi: Cone of influence
e-folding time dofmin: Minimum degrees of freedom (2 for Morlet)

## Details

The Fourier factor allows conversion between wavelet scale s and
equivalent Fourier period T: T = (4\*pi / (k_0 + sqrt(2 + k_0^2))) \* s

The cone of influence (COI) defines the e-folding time for edge effects:
COI = fourier_factor / sqrt(2)

Values outside the COI are subject to edge artifacts and should be
interpreted with caution.

## References

Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet
Analysis. Bulletin of the American Meteorological Society, 79(1), 61-78.
See Table 1 and Section 3f.

## Examples

``` r
# Standard Morlet wavelet parameters
params <- morlet_parameters(k0 = 6)
print(params)
#> fourier_factor            coi         dofmin 
#>      1.0330436      0.7304722      2.0000000 
# fourier_factor: 1.03
# coi: 0.73
# dofmin: 2
```
