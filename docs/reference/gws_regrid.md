# Regrid GWS to target period vector

Interpolates GWS from wavelet analysis onto a target period grid. Uses
linear interpolation and fills edge NAs.

## Usage

``` r
gws_regrid(wv, target_period, use_unmasked = FALSE)
```

## Arguments

- wv:

  List output from wavelet_spectral_analysis()

- target_period:

  Numeric vector of target periods

- use_unmasked:

  Logical. If TRUE and available, use gws_unmasked

## Value

Numeric vector of GWS values on target_period grid
