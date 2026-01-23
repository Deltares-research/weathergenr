# Additive wavelet analysis with CWT diagnostics

Combines CWT analysis (visualization and significance testing) with
MODWT MRA (additive components). CWT is used for diagnostics only; MODWT
MRA provides the additive component extraction used downstream.

This function expects
[`analyze_wavelet_spectrum()`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md)
to be available elsewhere in the package.

## Usage

``` r
analyze_wavelet_additive(
  series,
  signif = 0.9,
  noise = "red",
  min_period = 2,
  detrend = FALSE,
  filter = c("la8", "haar", "d4", "d6", "d8", "la16"),
  n_levels = NULL,
  boundary = "periodic",
  include_smooth = TRUE,
  max_period_frac = 1/3,
  cwt_mode = c("fast", "complete"),
  diagnostics = FALSE
)
```

## Arguments

- series:

  Numeric vector. Input series (regularly spaced, no missing values).

- signif:

  Numeric scalar in (0, 1). Significance level for CWT.

- noise:

  Character. Background noise model for CWT significance testing.
  Typically `"white"` or `"red"`.

- min_period:

  Numeric scalar. Minimum Fourier period for CWT.

- detrend:

  Logical. If TRUE, removes a linear trend before CWT.

- filter:

  Character. MODWT filter name. One of `"la8"` (default), `"haar"`,
  `"d4"`, `"d6"`, `"d8"`, `"la16"`.

- n_levels:

  Integer scalar or NULL. MODWT levels. If NULL, chosen from CWT
  significance (when available) and constrained by stability and
  max-period caps.

- boundary:

  Character. Boundary handling method. Only `"periodic"` is supported.

- include_smooth:

  Logical. If TRUE, includes the smooth component SJ.

- max_period_frac:

  Numeric scalar in (0, 1\]. No structure beyond this fraction of n.
  Default is `1/3`.

- cwt_mode:

  Character. CWT routine mode. One of `"fast"` or `"complete"`.

- diagnostics:

  Logical. If TRUE, attaches additional diagnostics (caps and
  covariance).

## Value

A list with class `"wavelet_additive"`. The object contains: - `cwt`:
output from
[`analyze_wavelet_spectrum()`](https://deltares-research.github.io/weathergenr/reference/analyze_wavelet_spectrum.md) -
`components`: additive MODWT MRA components (matrix) -
`component_names`: component column names - `periods`: representative
MODWT periods per component - `variance`, `variance_fraction`: component
variance and fraction of total - `n_levels`, `filter`, `filter_length`,
`boundary`, `include_smooth`, `max_period_frac` - `has_significance`,
`cwt_signif_periods`: CWT significance flags and periods -
`cwt_to_modwt_map`: MODWT levels implied by significant CWT periods
(under caps) - variance-accounting fields copied from
[`modwt_mra()`](https://deltares-research.github.io/weathergenr/reference/modwt_mra.md) -
optional `diagnostics` list when requested
