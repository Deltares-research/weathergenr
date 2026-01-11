# Extract significance curve from wavelet analysis result

Searches for GWS significance curve in wavelet analysis output. Tries
multiple possible field names.

## Usage

``` r
extract_signif_curve(wv)
```

## Arguments

- wv:

  List output from wavelet_spectral_analysis()

## Value

Numeric vector of significance values, or NULL if not found
