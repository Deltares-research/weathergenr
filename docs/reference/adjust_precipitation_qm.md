# Adjust Daily Precipitation with Gamma Quantile Mapping Under Monthly Scenario Factors

Applies a month-wise \*\*Gamma quantile mapping (QM)\*\* to daily
precipitation to impose prescribed changes in \*\*wet-day intensity\*\*
(monthly wet-day mean and variance).

The workflow is:

1.  Identify "wet-day intensities" as `precip > intensity_threshold`.

2.  For each calendar month, fit a \*\*baseline Gamma distribution\*\*
    to wet-day intensities using
    [`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).

3.  For each **simulation year index** `y` and month `m`, construct a
    \*\*target Gamma distribution\*\* by scaling baseline moments using
    `mean_factor[y, m]` and `var_factor[y, m]` (and optionally
    `scale_var_with_mean`).

4.  Map each wet-day intensity through the baseline CDF (`pgamma`) and
    then through the target inverse CDF (`qgamma`).

Values `<= intensity_threshold` are treated as "dry" for this function
and are returned unchanged. Therefore, this function does **not** change
wet-day frequency. If you need changes in wet/dry spell structure, apply
a separate occurrence perturbation.

## Usage

``` r
adjust_precipitation_qm(
  precip = NULL,
  mean_factor = NULL,
  var_factor = NULL,
  scale_var_with_mean = TRUE,
  exaggerate_extremes = FALSE,
  extreme_prob_threshold = 0.95,
  extreme_k = 1.2,
  enforce_target_mean = TRUE,
  month = NULL,
  year = NULL,
  intensity_threshold = 0,
  fit_method = "mme",
  min_events = 10,
  validate_output = TRUE,
  diagnostics = FALSE,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- precip:

  Numeric vector. Daily precipitation (mm/day). Must be `>= 0`. `NA`
  values are allowed and pass through unchanged.

- mean_factor:

  Numeric matrix with dimensions `(n_years x 12)`. Monthly
  multiplicative factors applied to the **baseline wet-day mean** for
  each month. Entry `mean_factor[y, m]` is used for simulation year
  index `y` and calendar month `m`.

- var_factor:

  Numeric matrix with dimensions `(n_years x 12)`. Monthly
  multiplicative factors applied to the **baseline wet-day variance**
  for each month. Entry `var_factor[y, m]` is used for simulation year
  index `y` and calendar month `m`.

- scale_var_with_mean:

  Logical scalar. If `TRUE`, the effective variance factor used to build
  target distributions is: \$\$var\\factor\\use = var\\factor \times
  mean\\factor^2\$\$ This tends to keep the coefficient of
  variation (CV) more stable when mean changes are applied. Default is
  `FALSE`.

- exaggerate_extremes:

  Logical scalar. If `TRUE`, applies a tail-exponent transform in
  probability space above `extreme_prob_threshold` before mapping into
  the target distribution. This amplifies relative changes in the upper
  tail while keeping the Gamma-to-Gamma mapping structure. Default is
  `FALSE`.

- extreme_prob_threshold:

  Numeric scalar in `(0, 1)`. Defines the start of the "extreme" tail in
  baseline probability space. For example, `0.95` corresponds to the top
  5% of wet-day intensities under the baseline month-specific Gamma.
  Used only when `exaggerate_extremes = TRUE`. Default is `0.95`.

- extreme_k:

  Numeric scalar `> 0`. Tail exponent controlling tail amplification
  when `exaggerate_extremes = TRUE`. Values `> 1` amplify extremes;
  values in `(0, 1)` dampen them. Default is `1.2`.

- enforce_target_mean:

  Logical scalar. If `TRUE`, rescales mapped wet-day values within each
  `(year index, month)` group so that the wet-day mean matches the
  intended target mean `baseline_mean(month) * mean_factor[y, m]`. This
  is most relevant when `exaggerate_extremes = TRUE` because tail
  amplification can shift the mean away from the intended target.
  Default is `TRUE`.

- month:

  Integer vector, same length as `precip`. Calendar month for each day
  (`1`–`12`).

- year:

  Integer vector, same length as `precip`. **Simulation year index** for
  each day (`1 =` first simulated year, `2 =` second, ...). Must be a
  contiguous index set `1:n_years`. Do not pass calendar years.
  `max(year)` must equal `nrow(mean_factor)` and `nrow(var_factor)`.

- intensity_threshold:

  Numeric scalar `>= 0`. Defines which values are treated as wet-day
  intensities for fitting and mapping: `precip > intensity_threshold`.
  Values `<= intensity_threshold` are returned unchanged. Default is
  `0`.

- fit_method:

  Character scalar. Estimation method passed to
  [`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
  for Gamma fitting (e.g., `"mme"`, `"mle"`). Default is `"mme"`.

- min_events:

  Integer scalar. Minimum number of wet-day intensities required within
  a month to fit the baseline Gamma. Months with fewer wet-day values
  are skipped; wet-day values in skipped months pass through unchanged.
  Default is `10`.

- validate_output:

  Logical scalar. If `TRUE`, replaces any `Inf`/`NaN` produced by the
  mapping with the original values at those positions and clamps any
  negative outputs to zero. Default is `TRUE`.

- diagnostics:

  Logical scalar. If `TRUE`, returns a list containing: the adjusted
  series, diagnostics from
  [`validate_quantile_mapping()`](https://deltares-research.github.io/weathergenr/reference/validate_quantile_mapping.md),
  and the fitted objects needed by higher-level workflows (baseline and
  target Gamma parameters). Default is `FALSE`.

- seed:

  Optional integer. If provided, sets a temporary RNG seed via
  [`set.seed()`](https://rdrr.io/r/base/Random.html). The previous RNG
  state is restored on exit. Default is `NULL`.

- verbose:

  Logical scalar. If `TRUE`, prints progress and warnings for skipped
  months and failed fits. Default is `FALSE`.

## Value

If `diagnostics = FALSE`, returns a numeric vector of the same length as
`precip`. The returned vector has attributes:

- `perturbed_months`: months (1–12) successfully fitted and perturbed

- `skipped_months`: months (1–12) skipped due to insufficient data or
  failed fits

- `n_failed_fits`: number of monthly Gamma fits that failed

If `diagnostics = TRUE`, returns a list:

- `adjusted`: numeric vector as above (with the same attributes)

- `diagnostics`: output of
  [`validate_quantile_mapping()`](https://deltares-research.github.io/weathergenr/reference/validate_quantile_mapping.md)

## Details

**Interpretation of mean/variance factors**

- Factors apply to wet-day intensities only
  (`precip > intensity_threshold`).

- They do not control wet-day frequency; apply a separate occurrence
  model if required.

## See also

[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html),
[`pgamma`](https://rdrr.io/r/stats/GammaDist.html),
[`qgamma`](https://rdrr.io/r/stats/GammaDist.html),
[`validate_quantile_mapping`](https://deltares-research.github.io/weathergenr/reference/validate_quantile_mapping.md),
[`diagnose_precip_qm`](https://deltares-research.github.io/weathergenr/reference/diagnose_precip_qm.md)
