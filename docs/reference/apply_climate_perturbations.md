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

- temp_transient:

  Logical. If TRUE, temperature deltas ramp linearly from zero to twice
  \`temp_delta\` across simulation years, preserving the same mean
  change.

- precip_transient:

  Logical. If TRUE, precipitation mean and variance factors ramp
  linearly across simulation years using the same transient logic.

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

All precipitation perturbations rely on a simulation-year index
\`year_idx = 1..n_years\`, derived internally from \`date\`. Calendar
years are not passed downstream. Any factor matrices supplied as
\`n_years x 12\` are indexed using this simulation-year convention.

## See also

\[adjust_precipitation_qm()\], \[calculate_monthly_pet()\]
