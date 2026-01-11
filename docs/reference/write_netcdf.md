# Write Gridded Data to a NetCDF File (Template-Based)

Writes gridded climate/weather time series to a NetCDF file using a
template NetCDF as the schema source (dimensions, spatial reference
variable, and global attributes). The function expects a list of
grid-cell payloads (each containing identical variable names and
equal-length time series) and places each cell into the output raster
cube according to `coord.grid` indices. Output is written as NetCDF-4
with optional compression and optional early rounding of values.

## Usage

``` r
write_netcdf(
  data = NULL,
  coord.grid = NULL,
  output.path = NULL,
  origin.date = NULL,
  calendar.type = "noleap",
  nc.template.file = NULL,
  nc.compression = 4,
  nc.spatial.ref = "spatial_ref",
  nc.file.prefix = "clim_change_rlz",
  nc.file.suffix = "",
  signif.digits = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  List of length `nrow(coord.grid)`. Each element is a named list of
  numeric vectors (time series) for each variable to be written. All
  list elements must have identical variable names and identical time
  series lengths.

- coord.grid:

  Data frame mapping grid cells to template indices. Must contain
  integer columns `xind` and `yind` (1-based indices) and have
  `nrow(coord.grid) == length(data)`. Indices must be within the
  template x/y bounds.

- output.path:

  Character string. Directory where the NetCDF file will be written.
  Created if it does not exist.

- origin.date:

  Date or character coercible to `Date`. Used to define the time units
  (days since origin date at 00:00:00).

- calendar.type:

  Character string. NetCDF calendar attribute for the time dimension
  (e.g., `"noleap"`). Passed to
  [`ncdf4::ncdim_def()`](https://rdrr.io/pkg/ncdf4/man/ncdim_def.html).

- nc.template.file:

  Character string. Path to an existing template NetCDF file used to
  infer dimensions, coordinate values, spatial reference, and global
  attributes.

- nc.compression:

  Integer in `[0, 9]`. NetCDF-4 deflate compression level applied to
  output variables.

- nc.spatial.ref:

  Character string. Name of the spatial reference variable in the
  template to copy to output (default `"spatial_ref"`). The variable
  must contain attributes `x_dim` and `y_dim`.

- nc.file.prefix:

  Character string. Output file prefix (base name without extension).

- nc.file.suffix:

  Character string. Optional suffix appended to the file name after an
  underscore. If `""`, no suffix is added.

- signif.digits:

  Optional integer. If provided, all values in `data` are rounded using
  `round(x, signif.digits)` prior to writing.

- verbose:

  Logical scalar. If `TRUE`, progress messages are emitted using
  [`logger::log_info()`](https://daroczig.github.io/logger/reference/log_level.html)
  with a `[NetCDF]` tag.

## Value

Invisibly returns the full path to the written NetCDF file (character
scalar).

## Details

The template file provided via `nc.template.file` is opened to:

- validate that the spatial reference variable `nc.spatial.ref` exists;

- read `x_dim` and `y_dim` attributes from the spatial reference
  variable to determine the x/y dimension names in the template;

- copy the spatial reference variable values and its attributes into the
  output file;

- copy template global attributes into the output file (best-effort);

- define output x/y dimension coordinate values based on template
  dimension values.

A `time` dimension is defined as `0:(nt - 1)` with units
`"days since <origin.date> 00:00:00"` and the specified `calendar.type`.
Each variable in `data[[1]]` is written as a 3D array with dimension
order `time x y x x`. Grid cells are inserted at `[ , yind, xind ]`
locations.

The output filename is constructed as:
`file.path(output.path, paste0(nc.file.prefix, optional_suffix, ".nc"))`
where `optional_suffix` is `paste0("_", nc.file.suffix)` when
`nc.file.suffix != ""`.

## File structure

- Dimensions: `time`, template `y`, template `x`

- Variables: one variable per name in `names(data[[1]])`, plus
  `nc.spatial.ref`

- Global attributes: copied from template (best-effort) plus
  `creation_date` and `created_with`

## Side effects

Creates directories (if needed), writes a NetCDF file, and may overwrite
an existing output file with the same name (a warning is issued).

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal example (schematic)
coord.grid <- data.frame(xind = c(1, 2), yind = c(1, 1))
data <- list(
  list(prcp = rnorm(365), tavg = rnorm(365)),
  list(prcp = rnorm(365), tavg = rnorm(365))
)

out_file <- write_netcdf(
  data = data,
  coord.grid = coord.grid,
  output.path = "out",
  origin.date = "1990-01-01",
  calendar.type = "noleap",
  nc.template.file = "template.nc",
  nc.file.prefix = "clim_change_rlz",
  nc.file.suffix = "run01",
  nc.compression = 4,
  signif.digits = 3,
  verbose = TRUE
)
} # }
```
