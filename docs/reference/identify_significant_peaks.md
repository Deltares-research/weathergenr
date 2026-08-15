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

## Details

Gaps in `gws` are interpolated; gaps in `gws_signif` are not, and the
asymmetry is deliberate. The spectrum is a measured quantity, so
interpolating across a gap in it is reasonable. The significance curve
is `NA` precisely where the cone of influence leaves too few effective
degrees of freedom for the test to run, and filling it there would
invent a test the record cannot support. A scale whose threshold is `NA`
yields no significant peak, which is the correct answer rather than a
missing one.

This matters most on short records. A 20-year annual series can support
a test only out to periods of roughly 3.5 years, so pass `gws_signif`
unfilled – supplying a filled curve will report peaks at decadal periods
from a threshold carried over from 3.5 years.
