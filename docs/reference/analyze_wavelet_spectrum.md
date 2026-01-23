# Continuous Wavelet Spectral Analysis with COI-Aware Significance

Performs continuous wavelet transform (CWT) analysis using a Morlet
mother wavelet (Torrence & Compo, 1998) and returns wavelet power, cone
of influence (COI), global wavelet spectrum (GWS), and significance
thresholds. The implementation includes:

- COI-aware global spectrum estimation,

- effective sample size (`neff`) based significance testing,

- AR(1) red-noise background estimation (Yule-Walker),

- optional scale-component reconstruction in `mode = "complete"`.

This function is intended for detection and diagnostic analysis of
dominant time scales. Reconstructed scale components are additive but
not orthogonal. For truly orthogonal decomposition suitable for ARIMA
modeling, see `analyze_wavelet_orthogonal` which combines CWT
diagnostics with MODWT-based orthogonal component extraction.

## Usage

``` r
analyze_wavelet_spectrum(
  series,
  signif = 0.9,
  noise = "red",
  min_period = 2,
  detrend = FALSE,
  mode = c("fast", "complete"),
  diagnostics = FALSE,
  lag1_ci = FALSE,
  lag1_ci_level = 0.95,
  lag1_boot_n = 500,
  seed = NULL,
  warn_neff = FALSE,
  neff_warn_min = 5,
  neff_warn_frac = 0.6
)
```

## Arguments

- series:

  Numeric vector. Input time series (regularly spaced, no missing
  values). Minimum length is 16 observations.

- signif:

  Numeric scalar in (0, 1). Significance level for wavelet power and GWS
  tests. Default is 0.90.

- noise:

  Character. Background noise model for significance testing. One of
  `"white"` or `"red"` (default). For red noise, an AR(1) model is
  estimated.

- min_period:

  Numeric scalar. Minimum Fourier period (in time units) used when
  identifying significant scales (default: 2).

- detrend:

  Logical. If `TRUE`, removes a linear-trend slope only (mean preserved)
  prior to analysis (default: `FALSE`).

- mode:

  Character. `"fast"` (default) returns essential outputs only;
  `"complete"` additionally returns pointwise significance and
  (optionally) reconstruction products.

- diagnostics:

  Logical. If `TRUE`, returns an additional `diagnostics` list. When
  `mode = "complete"`, diagnostics also include reconstruction error
  metrics when reconstructed components are computed.

- lag1_ci:

  Logical. If `TRUE`, computes a bootstrap confidence interval for the
  AR(1) coefficient (diagnostic only; does not affect inference).

- lag1_ci_level:

  Numeric scalar in (0, 1). Confidence level for the AR(1) bootstrap
  interval.

- lag1_boot_n:

  Integer. Number of bootstrap replicates for the AR(1) interval (\>=
  50).

- seed:

  Optional numeric scalar. Random seed used for bootstrap diagnostics.

- warn_neff:

  Logical. If `TRUE`, warns when effective sample size is small for most
  scales.

- neff_warn_min:

  Numeric scalar. Threshold below which neff is considered small.

- neff_warn_frac:

  Numeric scalar in (0, 1). Fraction of scales with
  `neff < neff_warn_min` required to trigger a warning.

## Value

A list. Always includes (at minimum):

- gws:

  Numeric vector. COI-masked global wavelet spectrum.

- gws_unmasked:

  Numeric vector. Unmasked global wavelet spectrum (plotting only).

- period:

  Numeric vector. Fourier periods corresponding to wavelet scales.

- gws_signif:

  Numeric vector. COI- and neff-corrected GWS significance threshold
  (inference).

- gws_signif_unmasked:

  Numeric vector. Unmasked GWS significance threshold (plotting only).

- has_significance:

  Logical. Whether any significant scales were detected.

- signif_periods:

  Integer vector. Indices of significant scales retained (one per
  contiguous band).

- coi:

  Numeric vector. Cone of influence in time units.

- power:

  Numeric matrix. Wavelet power spectrum (scales x time).

In `mode = "complete"`, additional fields may include:

- power_coi:

  Numeric matrix. COI-masked wavelet power.

- sigm:

  Numeric matrix. Pointwise significance ratio.

- sigm_coi:

  Numeric matrix. COI-masked pointwise significance ratio.

- power_signif_coi:

  Logical matrix. COI-masked pointwise significance mask.

- wave:

  Complex matrix. Wavelet coefficients (scales x time).

- comps:

  Numeric matrix. Reconstructed significant-scale components plus
  residual (if computed).

- comps_names:

  Character vector. Component names.

- gws_n_coi:

  Numeric vector. Number of time points inside COI per scale.

- neff:

  Numeric vector. Effective sample size per scale (masked).

- neff_unmasked:

  Numeric vector. Effective sample size per scale (unmasked).

- diagnostics:

  List. Returned only when `diagnostics = TRUE`.

## References

Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61-78.

## See also

[`morlet_wavelet`](https://deltares-research.github.io/weathergenr/reference/morlet_wavelet.md),
[`extract_wavelet_components`](https://deltares-research.github.io/weathergenr/reference/extract_wavelet_components.md),
`analyze_wavelet_orthogonal`,
[`plot_wavelet_power`](https://deltares-research.github.io/weathergenr/reference/plot_wavelet_power.md),
[`plot_wavelet_global_spectrum`](https://deltares-research.github.io/weathergenr/reference/plot_wavelet_global_spectrum.md)
