# Write Gridded Data to a NetCDF File (Template-Based)

Writes gridded climate/weather time series to a NetCDF file using a
template NetCDF as the schema source (dimensions, spatial reference
variable, and global attributes).

## Usage

``` r
write_netcdf(
  data = NULL,
  grid = NULL,
  out_dir = NULL,
  origin_date = NULL,
  calendar = "noleap",
  dates = NULL,
  template_path = NULL,
  compression = 4,
  spatial_ref = "spatial_ref",
  file_prefix = "clim_change_rlz",
  file_suffix = "",
  signif_digits = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  List of length `nrow(grid)`. Each element is a named list of numeric
  vectors (time series) for each variable to be written. All list
  elements must have identical variable names and equal length.

- grid:

  data.frame with at least `xind` and `yind` integer indices mapping
  each list element to an output grid cell.

- out_dir:

  Character. Output directory.

- origin_date:

  Character or Date. **The date of the first time step**, not an
  arbitrary reference epoch. The time coordinate is written as
  `0:(nt - 1)`, so the axis is only correct when `origin_date` is the
  series start – passing a conventional epoch such as `"1970-01-01"` for
  a series that begins in 2020 silently relocates it by 50 years. Pass
  `dates` to have this checked rather than assumed.

- calendar:

  Character. CF calendar for the time axis. One of `"noleap"`
  (`"365_day"`), `"standard"`, `"gregorian"`, `"proleptic_gregorian"`,
  or `"julian"`. Must match the calendar the data are actually on:
  [`generate_weather`](https://deltares-research.github.io/weathergenr/reference/generate_weather.md)
  output is 365-day, so `"noleap"` is correct for it, and mislabelling
  it shifts every date for any CF-aware reader.

- dates:

  Optional Date vector, one entry per time step, giving the axis the
  data are on. When supplied it is validated against `nt`,
  `origin_date`, and `calendar`, and the time coordinate is computed
  from it instead of being assumed contiguous. When `NULL` (default) the
  axis is assumed to be `nt` consecutive days starting at `origin_date`.

- template_path:

  Character. Path to a template NetCDF file.

- compression:

  Integer 0-9. NetCDF4 deflation level.

- spatial_ref:

  Character. Spatial reference variable name in template.

- file_prefix:

  Character. Prefix for output filename.

- file_suffix:

  Character. Optional suffix appended to filename.

- signif_digits:

  Integer. If not `NULL`, round values to this many significant digits.

- verbose:

  Logical. If TRUE, emit progress logs via
  [`.log()`](https://deltares-research.github.io/weathergenr/reference/dot-log.md).

## Value

Invisibly returns the written file path.
