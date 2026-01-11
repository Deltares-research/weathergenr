# Plot Wavelet Power Spectrum and Global Wavelet Spectrum

Creates a two-panel figure showing (i) the time-period wavelet power
spectrum (with cone of influence and significance contour) and (ii) the
global wavelet spectrum for a time series.

## Usage

``` r
plot_wavelet_spectra(
  variable,
  variable.year = NULL,
  period,
  POWER,
  GWS,
  GWS_signif,
  coi,
  sigm,
  variable.unit = "mm"
)
```

## Arguments

- variable:

  Numeric vector. The original time series (for length).

- variable.year:

  Numeric vector. Year labels for x-axis (optional).

- period:

  Numeric vector. Periods (scales) used in the wavelet transform.

- POWER:

  Matrix. Wavelet power spectrum (periods x time).

- GWS:

  Numeric vector. Global wavelet spectrum (mean power at each period).

- GWS_signif:

  Numeric vector. Significance threshold for the global wavelet
  spectrum.

- coi:

  Numeric vector. Cone of influence for each time point.

- sigm:

  Matrix. Significance mask or ratio (periods x time). Contour is drawn
  at 1.

- variable.unit:

  Character. Unit label for the variable (default is "mm").

## Value

A patchwork plot object with the time-period power spectrum and global
wavelet spectrum.

## Details

Key plotting choices for interpretability:

- Power is plotted on a log2 scale with robust color limits (5th-95th
  percentiles) to avoid a "washed out" field when the dynamic range is
  large.

- The spectrum is plotted on the native wavelet grid (no interpolation),
  avoiding smoothing/extrapolation artifacts that can flatten contrast.

- Uses a perceptually uniform `viridis` color scale.
