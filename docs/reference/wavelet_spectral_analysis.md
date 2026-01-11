# Continuous Wavelet Spectral Analysis with COI-Aware Significance

Performs continuous wavelet transform (CWT) analysis using the Morlet
mother wavelet following Torrence & Compo (1998), with extensions for:

- cone-of-influence (COI) aware global spectrum estimation,

- effective sample size (neff)-based significance testing,

- defensible AR(1) red-noise background estimation,

- stable and complete signal reconstruction from significant scales.

The function is designed for \*\*detection and diagnostic analysis of
dominant time scales\*\*, not for producing orthogonal components.
Reconstructed components are scale-localized but not independent.

## Usage

``` r
wavelet_spectral_analysis(
  variable,
  signif.level = 0.9,
  noise.type = "red",
  period.lower.limit = 2,
  detrend = FALSE,
  mode = c("fast", "complete"),
  return_diagnostics = TRUE,
  return_recon_error = FALSE,
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

- variable:

  Numeric vector. Input time series (regularly spaced, no missing
  values). Minimum length is 16 observations.

- signif.level:

  Numeric scalar in (0, 1). Significance level for wavelet power and
  global wavelet spectrum tests. Default is 0.90.

- noise.type:

  Character. Background noise model for significance testing. One of
  `"white"` or `"red"` (default). For red noise, an AR(1) model is
  estimated using Yule-Walker.

- period.lower.limit:

  Numeric scalar. Minimum Fourier period (in time units) to consider
  when identifying significant scales. Default is 2.

- detrend:

  Logical. If `TRUE`, removes a linear trend \*\*slope only\*\*
  (preserves the series mean level) before wavelet analysis. Default is
  `FALSE`.

- mode:

  Character. Computation mode: `"fast"` (default) for speed-optimized
  filtering with essential outputs only, or `"complete"` for
  comprehensive analysis including reconstruction and diagnostics. Fast
  mode returns 9 outputs (~50 KB), complete mode returns 18+ outputs (~5
  MB). Fast mode is ~3x faster by skipping signal reconstruction and
  detailed diagnostics.

- return_diagnostics:

  Logical. DEPRECATED. Use `mode = "complete"` instead. If `TRUE`,
  forces complete mode and returns internal diagnostic quantities.

- return_recon_error:

  Logical. If `TRUE`, returns scalar reconstruction error metrics
  verifying closure of reconstructed components.

- lag1_ci:

  Logical. If `TRUE`, computes a bootstrap confidence interval for the
  AR(1) coefficient (diagnostic only; does not affect inference).

- lag1_ci_level:

  Numeric scalar in (0, 1). Confidence level for the AR(1) bootstrap
  interval. Default is 0.95.

- lag1_boot_n:

  Integer. Number of bootstrap replicates used to estimate the AR(1)
  confidence interval. Default is 500.

- seed:

  Optional numeric scalar. Random seed for reproducible bootstrap
  diagnostics. Default is `NULL`.

- warn_neff:

  Logical. If `TRUE` (default), emits a warning when the effective
  sample size is small for most scales.

- neff_warn_min:

  Numeric scalar. Threshold below which neff is considered small.
  Default is 5.

- neff_warn_frac:

  Numeric scalar in (0, 1). Fraction of scales with
  `neff < neff_warn_min` required to trigger a warning. Default is 0.60.

## Value

A list with outputs depending on mode:

**Fast mode (9 outputs):**

- gws:

  Numeric vector. COI-masked global wavelet spectrum.

- gws_unmasked:

  Numeric vector. Unmasked global wavelet spectrum (plotting only).

- gws_period:

  Numeric vector. Fourier periods corresponding to wavelet scales.

- gws_signif:

  Numeric vector. COI- and neff-corrected significance threshold
  (inference).

- gws_signif_unmasked:

  Numeric vector. Unmasked significance threshold (plotting only).

- has_significance:

  Logical. Indicates whether any significant scales were found.

- signif_periods:

  Integer vector. Indices of significant wavelet scales.

- coi:

  Numeric vector. Cone of influence in time units.

- power:

  Numeric matrix. Wavelet power spectrum (for basic plotting).

**Complete mode (additional 10+ outputs):**

- power_coi:

  Numeric matrix. COI-masked wavelet power.

- sigm:

  Numeric matrix. Pointwise significance ratio (unmasked).

- sigm_coi:

  Numeric matrix. COI-masked pointwise significance ratio.

- power_signif_coi:

  Logical matrix. COI-masked pointwise significance mask.

- wave:

  Complex matrix. Wavelet coefficients (scales x time).

- comps:

  Numeric matrix. Reconstructed components at significant scales plus
  residual.

- comps_names:

  Character vector. Column names of comps.

- gws_n_coi:

  Numeric vector. Number of time points inside the COI per scale.

- gws_neff:

  Numeric vector. Effective sample size per scale (masked).

- gws_neff_unmasked:

  Numeric vector. Effective sample size per scale (unmasked).

- diagnostics:

  List. Returned only if `return_diagnostics = TRUE`.

- reconstruction_error:

  List. Returned only if `return_recon_error = TRUE`.

## Details

The implementation follows Torrence & Compo (1998) for the Morlet
wavelet, with the following methodological refinements:

- two-sided (TC98-consistent) wavelet normalization,

- COI-masked global wavelet spectrum,

- effective degrees of freedom based on wavelet decorrelation time,

- separation of inference curves from plotting curves.

Reconstructed components correspond to \*\*individual wavelet scales\*\*
and are additive but not orthogonal. They should not be interpreted as
independent stochastic modes.

## References

Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61-78.

## See also

[`morlet_wavelet`](https://deltares-research.github.io/weathergenr/reference/morlet_wavelet.md),
[`extract_wavelet_components`](https://deltares-research.github.io/weathergenr/reference/extract_wavelet_components.md)
