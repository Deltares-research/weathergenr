# Extract Global-Spectrum Significance Curve

Extracts a global-spectrum significance curve from a wavelet analysis
output list. Intended as a small compatibility helper for downstream
code that expects a single numeric vector. If not found, returns `NULL`.

## Usage

``` r
extract_signif_curve(wavelet)
```

## Arguments

- wavelet:

  List output from
  [`analyze_wavelet_spectrum`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md).

## Value

Numeric vector of significance values, or `NULL` if not found.
