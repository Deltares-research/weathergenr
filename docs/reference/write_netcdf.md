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

  Character or Date. Origin date used for NetCDF time units (e.g.,
  `"1970-01-01"`).

- calendar:

  Character. NetCDF calendar type (e.g., `"noleap"`).

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
