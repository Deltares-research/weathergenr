# Plot Wavelet Power and Global Wavelet Spectrum

Produces a two-panel diagnostic figure for a 1D time series wavelet
analysis: (1) the time-period wavelet power field and (2) the global
wavelet spectrum (time-averaged power by period) with its significance
threshold.

## Usage

``` r
plot_wavelet_power(
  series,
  time = NULL,
  period,
  power,
  gws,
  gws_signif,
  coi,
  signif_mask,
  unit = "mm"
)
```

## Arguments

- series:

  Numeric vector. The time series used to define the time dimension
  (length must match the number of columns in `power`).

- time:

  Optional numeric vector of length `length(series)` used for the
  x-axis. Must be strictly increasing. If `NULL`, uses
  `seq_len(length(series))`.

- period:

  Numeric vector of wavelet periods (scales), typically in years (or
  days), strictly increasing.

- power:

  Numeric matrix of wavelet power with dimensions
  `length(period) x length(series)`.

- gws:

  Numeric vector. Global wavelet spectrum (mean power for each
  `period`). Must have length `length(period)`.

- gws_signif:

  Numeric vector. Significance threshold for `gws`. Must have length
  `length(period)`.

- coi:

  Numeric vector of length `length(series)` giving the cone of influence
  in the same units as `period`. Values outside the plotted period range
  are ignored.

- signif_mask:

  Numeric matrix with dimensions `length(period) x length(series)`. A
  contour is drawn at `1`. This is typically a significance ratio or a
  binary mask encoded as 0/1 (or NA).

- unit:

  Character. Unit label used in the global spectrum axis label (default:
  `"mm"`).

## Value

A `patchwork` object combining the power-field panel and the global
spectrum panel.

## Details

**What this plot is for**

- Identify when specific periodicities are active in time (power field).

- Summarize dominant variability scales across the full record (global
  spectrum).

- Check edge effects and statistical significance (COI and contour).

**Key plotting choices**

- Power is shown on a \\\log_2\\ scale after discarding non-finite
  values and non-positive power.

- The power field is drawn on the native wavelet grid (no interpolation)
  using explicit cell boundaries, which avoids smoothing artifacts.

- Color limits are set using robust quantiles (5th-95th percentile) to
  avoid low-contrast fields when dynamic range is large.

**Inputs must be consistent**

- `n_time = length(series)` and `n_period = length(period)`.

- `power` and `signif_mask` must have dimensions `n_period x n_time`.

- `coi` must have length `n_time`.

- `period` must be strictly increasing.

- If supplied, `time` must be strictly increasing and have length
  `n_time`.

## See also

[`plot_wavelet_global_spectrum`](https://deltares-research.github.io/weathergenr/reference/plot_wavelet_global_spectrum.md)
for a standalone global spectrum plot.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you computed wavelet outputs elsewhere:
# w <- analyze_wavelet_spectrum(series)
p <- plot_wavelet_power(
  series      = series,
  time        = years,
  period      = w$period,
  power       = w$power,
  gws         = w$gws,
  gws_signif  = w$gws_signif,
  coi         = w$coi,
  signif_mask = w$signif_mask,
  unit        = "mm"
)
print(p)
} # }
```
