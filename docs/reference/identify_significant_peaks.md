# Identify significant peaks in an observed global wavelet spectrum

Detects local maxima in an observed global wavelet spectrum and keeps
only those peaks whose power exceeds the corresponding significance
curve value.

Peaks are ranked by signal to noise ratio defined as power divided by
significance, then by power. Up to `n_max` peaks are returned.

If the significance curve is missing or is not aligned to `gws`, the
function returns an empty result, treating the significance information
as unavailable.

## Usage

``` r
identify_significant_peaks(gws, gws_signif, period, n_max = 3L)
```

## Arguments

- gws:

  Numeric vector. Observed global wavelet spectrum values.

- gws_signif:

  Numeric vector or NULL. Significance curve aligned to `gws`. Must have
  the same length as `gws`.

- period:

  Numeric vector. Period values associated with `gws`. Must have the
  same length as `gws`.

- n_max:

  Integer scalar. Maximum number of significant peaks to return.

## Value

Data frame with one row per selected peak and columns: `idx` index in
the spectrum, `period` the period at the peak, `power` the peak power,
`signif` the significance curve value at the peak, `snr` signal to noise
ratio. Returns a data frame with zero rows if no significant peaks are
found.
