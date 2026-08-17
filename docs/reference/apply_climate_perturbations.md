# Apply Monthly Climate Perturbations to Gridded Daily Weather Series

Applies monthly climate perturbations to daily gridded weather
simulations in a modular and physically interpretable way. The function
operates independently on precipitation intensity, precipitation
occurrence, temperature, and (optionally) potential evapotranspiration
(PET).

The perturbation workflow is: 1. Construct a simulation-year index from
\`date\` (\`1..n_years\`). 2. Optionally modify monthly wet-day
probabilities within each (year index, month) group. 3. Apply
precipitation intensity perturbations to wet days using Gamma-based
quantile mapping via \[adjust_precipitation_qm()\]. 4. Apply temperature
perturbations using monthly additive deltas (step change or transient
ramp). 5. Optionally recompute PET from perturbed temperature fields.

This function is intended for climate stress testing and scenario
analysis, where controlled changes in mean climate, variability, and
extremes are required while preserving realistic daily structure.

## Usage

``` r
apply_climate_perturbations(
  data = NULL,
  grid = NULL,
  date = NULL,
  precip_mean_factor = NULL,
  precip_var_factor = NULL,
  precip_occurrence_factor = NULL,
  precip_intensity_threshold = 0,
  temp_delta = NULL,
  temp_range_factor = NULL,
  temp_transient = TRUE,
  precip_transient = TRUE,
  precip_occurrence_transient = TRUE,
  compute_pet = TRUE,
  pet_method = "hargreaves",
  qm_fit_method = "mme",
  scale_var_with_mean = TRUE,
  exaggerate_extremes = FALSE,
  extreme_prob_threshold = 0.95,
  extreme_k = 1.2,
  enforce_target_mean = TRUE,
  precip_cap_mm_day = NULL,
  precip_floor_mm_day = NULL,
  precip_cap_quantile = NULL,
  seed = NULL,
  verbose = FALSE,
  diagnostic = TRUE
)
```

## Arguments

- data:

  List of data.frames, one per grid cell. Each data.frame must contain
  \`precip\`, \`temp\`, \`temp_min\`, and \`temp_max\`. If \`compute_pet
  = TRUE\`, a column \`pet\` is added or overwritten.

- grid:

  Data frame of grid metadata with \`nrow(grid) == length(data)\`. Must
  include a latitude column: \`lat\` (preferred) or \`y\` (legacy).

- date:

  Date vector with length equal to the number of rows in each grid-cell
  data.frame. Used to derive month and simulation-year indices.

- precip_mean_factor:

  Numeric vector of length 12 or numeric matrix \`(n_years x 12)\`.
  Monthly multiplicative factors applied to wet-day precipitation mean
  during intensity perturbation.

- precip_var_factor:

  Numeric vector of length 12 or numeric matrix \`(n_years x 12)\`.
  Monthly multiplicative factors applied to wet-day precipitation
  variance.

- precip_occurrence_factor:

  Optional numeric vector of length 12 or numeric matrix \`(n_years x
  12)\`. Multiplicative factors applied to monthly wet-day probability.

- precip_intensity_threshold:

  Numeric scalar \`\>= 0\`. Defines wet days as \`precip \>
  precip_intensity_threshold\`. Also passed to
  \[adjust_precipitation_qm()\] as \`intensity_threshold\`. Default is
  0.

- temp_delta:

  Numeric vector of length 12. Monthly additive temperature deltas
  (degC) applied to \`temp\`, \`temp_min\`, and \`temp_max\`.

- temp_range_factor:

  Optional numeric vector of length 12 or numeric matrix \`(n_years x
  12)\`. Multiplicative factors applied to the diurnal temperature
  range, \`temp_max - temp_min\`, scaled about its own midpoint so the
  daily mean is unchanged. \`NULL\` (default) leaves the range as it is.

  Worth setting deliberately rather than leaving alone. \`temp_delta\`
  is added to \`temp\`, \`temp_min\` and \`temp_max\` alike, so warming
  on its own leaves the diurnal range untouched — and the Hargreaves PET
  used here goes as the square root of that range. Warming alone
  therefore moves PET by only about 2.3 percent per degree, while the
  range is the term that decides how much more: at +4 degC, PET rises
  3.7 percent with the range a tenth smaller, 9.3 percent with it fixed,
  and 14.7 percent with it a tenth larger. Leaving this \`NULL\` is a
  choice for the middle of that span, not a neutral default, and it
  matters wherever PET drives the water balance.

  Observed and projected diurnal-range trends vary by region and season,
  so there is no defensible non-null default; supply values for the
  region being studied, or vary the factor as one of the stress-test
  dimensions.

- temp_transient:

  Logical. If TRUE, temperature deltas ramp linearly from zero to twice
  \`temp_delta\` across simulation years, preserving the same mean
  change. Also governs \`temp_range_factor\`, which ramps from 1 on the
  same principle.

- precip_transient:

  Logical. If TRUE, precipitation mean and variance factors ramp
  linearly across simulation years using the same transient logic: from
  1 to \`2f - 1\`, so the ramp averages to \`f\` over the period. See
  the note below on what may be supplied alongside it.

- precip_occurrence_transient:

  Logical. If TRUE, precipitation occurrence factors ramp linearly
  across simulation years.

- compute_pet:

  Logical. If TRUE, recompute PET from perturbed temperature using
  \[calculate_monthly_pet()\].

- pet_method:

  Character string specifying the PET method passed to
  \[calculate_monthly_pet()\] (default: "hargreaves").

- qm_fit_method:

  Character string specifying the Gamma fitting method passed to
  \[adjust_precipitation_qm()\] (e.g., "mme", "mle").

- scale_var_with_mean:

  Logical. If TRUE, precipitation variance scaling is computed as
  \`var_factor_use = precip_var_factor \* precip_mean_factor^2\`.

- exaggerate_extremes:

  Logical. If TRUE, amplify upper-tail wet-day precipitation during the
  intensity quantile mapping step (forwarded to
  \[adjust_precipitation_qm()\]).

- extreme_prob_threshold:

  Numeric scalar in (0, 1). Probability threshold defining the start of
  the tail region for amplification in \[adjust_precipitation_qm()\].

- extreme_k:

  Numeric scalar \> 0. Tail exponent controlling amplification strength
  in \[adjust_precipitation_qm()\]. Values \> 1 amplify extremes; values
  in (0, 1) dampen.

- enforce_target_mean:

  Logical. If TRUE, rescale mapped wet-day values within each (year
  index, month) group so the wet-day mean matches the intended target
  mean after tail amplification.

- precip_cap_mm_day:

  Optional numeric scalar. Absolute upper cap (mm/day) applied to
  precipitation after all perturbations.

- precip_floor_mm_day:

  Optional numeric scalar. Minimum precipitation amount (mm/day) applied
  to wet days after perturbation.

- precip_cap_quantile:

  Optional numeric scalar in (0, 1). Quantile-based cap computed from
  the perturbed precipitation series and applied as an additional upper
  bound.

- seed:

  Optional integer. Sets the random seed for reproducible precipitation
  occurrence perturbations and any stochastic components.

- verbose:

  Logical. If TRUE, prints progress messages and warnings.

- diagnostic:

  Logical. If TRUE (default), return precipitation quantile-mapping
  diagnostics from \[adjust_precipitation_qm()\] alongside the perturbed
  data.

## Value

If \`diagnostic = TRUE\` (default), a list with: \* \`data\`: list of
data.frames (one per grid cell) \* \`diagnostic\`: list of per-grid
diagnostic objects returned by \[adjust_precipitation_qm()\], with the
large \`adjusted\` vector removed

If \`diagnostic = FALSE\`, only the list of data.frames is returned.

## Year indexing convention (critical)

All precipitation perturbations rely on a year index \`year_idx =
1..n_years\`, derived internally from \`date\` as \`calendar_year -
min(calendar_year) + 1\`. The absolute calendar year is not passed
downstream, but the \*counting\* is by calendar year, so \`n_years\`
here means the number of distinct calendar years spanned by \`date\` –
not the generator's \`n_years\` argument.

The two coincide only for a calendar-year run. Under a water year
(\`year_start_month \> 1\`), a 20-water-year series spans 21 calendar
years, and a factor matrix must then be \`21 x 12\`; a \`20 x 12\`
matrix is rejected. When in doubt, size matrices from the date vector
itself: \`length(unique(format(date, "%Y")))\`. Supplying a length-12
vector instead of a matrix sidesteps the question entirely, since it
applies to every year.

## Transient factors and year-varying matrices

A transient factor is specified by its \*\*end state\*\*, not by a
per-year path. \`transient = TRUE\` ramps linearly from 1 to \`2f - 1\`,
so the ramp averages to \`f\` over the period, and only the first row of
the factor is read – in transient mode one number per month is the
entire specification.

Supplying a year-varying \`n_years x 12\` matrix \*and\* setting the
matching transient flag is therefore a contradiction, and it is now an
error. Rows \`2:n_years\` would otherwise be discarded silently, which
is how a matrix encoding a 50 percent reduction could come back as no
change at all. Either pass a length-12 vector and let the ramp build the
path, or set the transient flag to \`FALSE\` and let the matrix be
applied year by year. A matrix whose rows are all identical is
equivalent to a vector and remains accepted.

## See also

\[adjust_precipitation_qm()\], \[calculate_monthly_pet()\]
