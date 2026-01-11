# Apply Climate Change Perturbations to Gridded Weather Data

Applies climate change perturbations to gridded daily weather data by
modifying precipitation and temperature using monthly change factors.
Perturbations can be applied as step changes or gradual transient
changes. Both approaches are constructed to yield the same mean change
over the simulation period.

## Usage

``` r
apply_climate_perturbations(
  climate.data = NULL,
  climate.grid = NULL,
  sim.dates = NULL,
  change.factor.precip.mean = NULL,
  change.factor.precip.variance = NULL,
  change.factor.temp.mean = NULL,
  transient.temp.change = TRUE,
  transient.precip.change = TRUE,
  calculate.pet = TRUE,
  fit.method = "mme",
  verbose = FALSE
)
```

## Arguments

- climate.data:

  List of data frames, one per grid cell. Each data frame must contain
  columns `precip`, `temp`, `temp_min`, `temp_max`, and optionally
  `pet`.

- climate.grid:

  Data frame containing grid cell metadata. Must include a column `y`
  representing latitude in decimal degrees.

- sim.dates:

  Vector of `Date` objects corresponding to the rows of each element in
  `climate.data`.

- change.factor.precip.mean:

  Numeric vector of length 12. Monthly multiplicative factors for
  precipitation mean representing the desired average change over the
  simulation period.

- change.factor.precip.variance:

  Numeric vector of length 12. Monthly multiplicative factors applied to
  precipitation variance.

- change.factor.temp.mean:

  Numeric vector of length 12. Monthly additive temperature changes in
  degrees Celsius representing the desired average change.

- transient.temp.change:

  Logical. If TRUE, temperature changes increase linearly from zero to
  twice the specified factor so that the target mean change is achieved.
  If FALSE, the full factor is applied uniformly.

- transient.precip.change:

  Logical. If TRUE, precipitation change factors increase linearly from
  1.0 to a calculated endpoint that yields the target mean. If FALSE,
  the full factor is applied uniformly.

- calculate.pet:

  Logical. If TRUE, potential evapotranspiration is recalculated from
  perturbed temperatures using the Hargreaves method.

- fit.method:

  Character string specifying the distribution fitting method passed to
  `quantile_mapping`.

- verbose:

  Logical. If TRUE, progress messages are printed.

## Value

A list of data frames with the same structure as `climate.data`,
containing perturbed weather variables. Invalid values are replaced with
zero.

## Details

Climate change perturbations are applied using the following mechanisms:

- Precipitation is adjusted using quantile mapping with mean and
  variance scaling.

- Temperature variables are adjusted using additive monthly shifts.

- Potential evapotranspiration is recalculated using the Hargreaves
  equation when requested.

Step and transient perturbations are formulated to produce the same mean
change over the full simulation period. Step changes apply the specified
factor uniformly across all years. Transient changes increase linearly
from a baseline to an endpoint that yields the same mean effect.

For precipitation, transient factors ramp from 1.0 to a computed upper
value. For temperature, transient changes ramp from zero to twice the
specified mean change. This ensures comparability between step and
transient scenarios while representing different temporal pathways.

During transient precipitation adjustments, a minimum factor of 0.01 is
enforced to prevent numerical instabilities.

## See also

[`quantile_mapping`](https://deltares-research.github.io/weathergenr/reference/quantile_mapping.md)

## Examples

``` r
if (FALSE) { # \dontrun{
perturbed_step <- apply_climate_perturbations(
  climate.data = weather_grids,
  climate.grid = grid_meta,
  sim.dates = dates,
  change.factor.precip.mean = rep(1.2, 12),
  change.factor.precip.variance = rep(1.0, 12),
  change.factor.temp.mean = rep(2.0, 12),
  transient.temp.change = FALSE,
  transient.precip.change = FALSE
)
} # }
```
