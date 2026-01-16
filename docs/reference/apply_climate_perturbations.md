# Apply monthly climate perturbations to gridded daily weather series

Applies \*\*monthly\*\* climate perturbations to \*\*daily\*\* gridded
weather series:

- \*\*Precipitation\*\* is perturbed using quantile mapping with monthly
  mean and variance factors via
  [`perturb_prcp_qm()`](https://deltares-research.github.io/weathergenr/reference/perturb_prcp_qm.md).

- \*\*Temperature\*\* (`temp`, `temp_min`, `temp_max`) is perturbed
  using monthly additive deltas (step or transient).

- \*\*Potential evapotranspiration (PET)\*\* can be recomputed from
  perturbed temperatures using
  [`calculate_monthly_pet()`](https://deltares-research.github.io/weathergenr/reference/calculate_monthly_pet.md).

Perturbations can be applied as a \*\*step change\*\* (constant in time)
or as a \*\*transient change\*\* (linearly ramping over years) while
preserving the same \*mean\* change over the simulation period.

## Usage

``` r
apply_climate_perturbations(
  data = NULL,
  grid = NULL,
  date = NULL,
  prcp_mean_factor = NULL,
  prcp_var_factor = NULL,
  temp_delta = NULL,
  temp_transient = TRUE,
  prcp_transient = TRUE,
  compute_pet = TRUE,
  pet_method = "hargreaves",
  qm_fit_method = "mme",
  verbose = FALSE
)
```

## Arguments

- data:

  List of data.frames, one per grid cell. Each data.frame must contain:

  - `prcp` (daily precipitation)

  - `temp` (daily mean temperature)

  - `temp_min` (daily minimum temperature)

  - `temp_max` (daily maximum temperature)

  If `compute_pet = TRUE`, `pet` is added or overwritten.

- grid:

  data.frame of grid metadata with `nrow(grid) == length(data)`. Must
  include column `lat` (latitude in decimal degrees).

- date:

  Date vector of length `nrow(data[[i]])` (identical across cells).

- prcp_mean_factor:

  Numeric vector length 12. Monthly multiplicative factors for
  precipitation mean.

- prcp_var_factor:

  Numeric vector length 12. Monthly multiplicative factors for
  precipitation variance.

- temp_delta:

  Numeric vector length 12. Monthly additive temperature deltas (degC)
  applied to `temp`, `temp_min`, and `temp_max`.

- temp_transient:

  Logical. If `TRUE`, temperature deltas ramp linearly from 0 to
  `2 * temp_delta` over years (mean equals `temp_delta`).

- prcp_transient:

  Logical. If `TRUE`, precipitation factors ramp linearly from 1 to
  `(factor - 1) * 2 + 1` over years.

- compute_pet:

  Logical. If `TRUE`, recompute PET using
  [`calculate_monthly_pet`](https://deltares-research.github.io/weathergenr/reference/calculate_monthly_pet.md).

- pet_method:

  Character. PET method passed to
  [`calculate_monthly_pet`](https://deltares-research.github.io/weathergenr/reference/calculate_monthly_pet.md)
  (default: `"hargreaves"`).

- qm_fit_method:

  Character. Distribution-fitting method for
  [`perturb_prcp_qm()`](https://deltares-research.github.io/weathergenr/reference/perturb_prcp_qm.md).

- verbose:

  Logical. Emit progress logs.

## Value

List of data.frames with perturbed variables.

## See also

[`perturb_prcp_qm`](https://deltares-research.github.io/weathergenr/reference/perturb_prcp_qm.md),
[`diagnose_prcp_qm`](https://deltares-research.github.io/weathergenr/reference/diagnose_prcp_qm.md),
[`calculate_monthly_pet`](https://deltares-research.github.io/weathergenr/reference/calculate_monthly_pet.md)
