# Read a gridded NetCDF file into tidy weather format (fast + memory efficient)

Efficiently reads 3D (x, y, time) gridded variables from a NetCDF file
and returns:

- per-grid-cell tidy time series (list of tibbles),

- grid coordinates and indices,

- a Date vector,

- dimension metadata,

- variable and global attributes.

Variable behaviour:

- `variables = NULL` -\> load all data variables (except `spatial.ref`).

- `variables = c("var1","var2")` -\> load only these variables.

- `var_rename = c(old1 = "new1", old2 = "new2")` -\> rename a subset of
  the loaded variables in the output (columns, attributes, dimnames).
  All names in `var_rename` must be among the selected variables.

## Usage

``` r
read_netcdf(
  nc.file,
  variables = NULL,
  var_rename = NULL,
  leap.days = TRUE,
  omit.empty = TRUE,
  spatial.ref = "spatial_ref",
  signif.digits = NULL
)
```

## Arguments

- nc.file:

  Character. Path to the NetCDF file.

- variables:

  Character vector of variable names to load. If `NULL`, all variables
  except `spatial.ref` are loaded. Default is `NULL`.

- var_rename:

  Named character vector giving new names for a subset of the loaded
  variables, e.g. `c(precip = "prcp", temp = "T")`. Old names must be a
  subset of `variables` (or of all variables if `variables=NULL`).
  Default is `NULL`.

- leap.days:

  Logical. Whether the data include leap days. If `FALSE`, Feb 29 is
  removed from the date vector. Default is `TRUE`.

- omit.empty:

  Logical. If `TRUE`, drop grid cells where all values are NA for all
  variables and all time steps. Default is `TRUE`.

- spatial.ref:

  Character. Name of the spatial reference variable in the NetCDF file
  (used to get x/y dimension names). Default is `"spatial_ref"`.

- signif.digits:

  Integer. If not `NULL`, values are rounded to this number of
  significant digits. Default is `NULL`.

## Value

A list with five components:

- data:

  List of tibbles, one per grid cell. Each tibble has rows representing
  time steps and columns representing variables.

- grid:

  Tibble with columns `id, xind, yind, x, y` describing the spatial grid
  cells.

- date:

  Vector of `Date` objects with length equal to the number of time
  steps.

- dimensions:

  Named list of dimension vectors (e.g., x, y, time) as extracted from
  the NetCDF file.

- attributes:

  Named list containing variable attributes (one element per variable
  using final renamed names), spatial reference attributes, and global
  attributes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load all variables from a NetCDF file
weather_data <- read_netcdf("climate_data.nc")

# Load only specific variables and rename them
weather_data <- read_netcdf(
  nc.file = "climate_data.nc",
  variables = c("precipitation", "temperature"),
  var_rename = c(precipitation = "prcp", temperature = "temp")
)

# Load without leap days
weather_data <- read_netcdf("climate_data.nc", leap.days = FALSE)

# Round values and keep only non-empty cells
weather_data <- read_netcdf(
  nc.file = "climate_data.nc",
  signif.digits = 4,
  omit.empty = TRUE
)
} # }
```
