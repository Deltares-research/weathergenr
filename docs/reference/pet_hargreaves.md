# Estimate Monthly Potential Evapotranspiration (Hargreaves Method)

Computes monthly mean potential evapotranspiration (PET) using the
Hargreaves temperature-based method. PET is estimated from mean air
temperature, diurnal temperature range, and extraterrestrial radiation
derived from latitude and time of year.

The implementation follows the FAO-56 formulation and returns PET in
units of millimeters per day (mm/day).

## Usage

``` r
pet_hargreaves(months, temp, tdiff, lat)
```

## Arguments

- months:

  Integer vector (1-12). Calendar months corresponding to each PET
  estimate.

- temp:

  Numeric vector. Mean monthly air temperature (degrees C).

- tdiff:

  Numeric vector. Monthly diurnal temperature range (Tmax ??? Tmin,
  degrees C). Must be non-negative.

- lat:

  Numeric scalar. Latitude in decimal degrees (positive north, negative
  south).

## Value

Numeric vector of monthly potential evapotranspiration values in mm/day,
with length equal to \`length(months)\`.

## Details

Extraterrestrial radiation is computed using a fixed representative day
of year for the middle of each month. Radiation is converted from MJ
m^2/day to equivalent evaporation depth using the FAO conversion factor
of 0.408.

The Hargreaves equation is: \$\$ PET = 0.0023 \* Ra \* (T + 17.8) \*
sqrt(Tmax - Tmin) \$\$ where Ra is extraterrestrial radiation.

## References

Hargreaves, G.H. and Samani, Z.A. (1985). Reference crop
evapotranspiration from temperature. Applied Engineering in Agriculture,
1(2), 96-99.

FAO (1998). Crop Evapotranspiration (FAO Irrigation and Drainage Paper
56).
