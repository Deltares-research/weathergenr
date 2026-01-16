# Plot Global Wavelet Spectrum

Plots the global wavelet spectrum (time-averaged wavelet power by
period) for an observed series, together with a period-wise significance
threshold. Optionally overlays the ensemble mean of simulated global
spectra.

## Usage

``` r
plot_wavelet_global_spectrum(period, signif, obs_power, sim_power = NULL)
```

## Arguments

- period:

  Numeric vector of wavelet periods (scales), typically in years (or
  days). Should be strictly increasing for meaningful interpretation.

- signif:

  Numeric vector of length `length(period)` giving the global-spectrum
  significance threshold at each period.

- obs_power:

  Numeric vector of length `length(period)` giving the observed global
  wavelet spectrum.

- sim_power:

  Optional numeric matrix of simulated global spectra with
  `nrow(sim_power) == length(period)` and one column per simulation.

## Value

A `ggplot` object.

## Details

**What this plot is for**

- Identify dominant variability scales (peaks in global power).

- Compare observed spectral structure to a simulated ensemble.

- Interpret variability relative to a significance threshold.

**Ensemble behavior** When `sim_power` is provided, the function adds
the ensemble mean as a solid line. This function intentionally keeps the
default view simple; if you want an ensemble envelope (e.g., min/max or
quantiles) add it outside this function.

## See also

[`plot_wavelet_power`](https://deltares-research.github.io/weathergenr/reference/plot_wavelet_power.md)
for a combined diagnostic plot including the time-period power field.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_wavelet_global_spectrum(
  period   = w$period,
  signif   = w$gws_signif,
  obs_power = w$gws,
  sim_power = w$gws_sim  # matrix: period x n_sim (optional)
)
print(p)
} # }
```
