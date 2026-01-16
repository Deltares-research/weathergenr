# Calculate Monthly Potential Evapotranspiration (PET)

Calculates monthly potential evapotranspiration (PET) using a selected
method. Currently supports the Hargreaves temperature-based method.

## Usage

``` r
calculate_monthly_pet(month, temp, temp_range, lat_deg, method = "hargreaves")
```

## Arguments

- month:

  Integer vector (1-12). Calendar month for each PET estimate.

- temp:

  Numeric vector. Mean monthly air temperature (degC).

- temp_range:

  Numeric vector. Monthly diurnal temperature range (Tmax - Tmin, degC).
  Must be non-negative.

- lat_deg:

  Numeric scalar. Latitude in decimal degrees (positive north).

- method:

  Character scalar. PET method. Currently `"hargreaves"`.

## Value

Numeric vector of monthly PET values in mm/day (length `length(month)`).
