#' Write Gridded Data to a NetCDF File (Template-Based)
#'
#' @description
#' Writes simulated or processed gridded climate/weather data to a NetCDF file, using variable and dimension
#' definitions from a template NetCDF. This wrapper ensures consistency with your reference files (dimensions, spatial ref, etc).
#'
#' Each element in `data` should be a list of variables for a single grid cell; the list should have the same structure and names as variables in the template file.
#'
#' @param data        A list of length = n_grids. Each element is a named list of vectors (or matrices) containing simulated data for each variable at that grid.
#'                   For example, \code{data[[i]]$precip} is a vector for grid i.
#' @param coord.grid  Data frame mapping grid indices (`xind`, `yind`) for each grid cell. Should have same length/order as `data`.
#' @param output.path Directory to write the file. Will be created if it doesn't exist.
#' @param origin.date Origin date as "YYYY-MM-DD" string for the time dimension.
#' @param nc.template.file Path to a reference/template NetCDF file (must exist).
#' @param calendar.type Calendar type for NetCDF time (default: "noleap").
#' @param nc.compression Integer; NetCDF compression level (0-9, default 4).
#' @param nc.spatial.ref Name of the spatial reference variable in the template (default "spatial_ref").
#' @param nc.file.prefix Prefix for the output file name (default "clim_change_rlz").
#' @param nc.file.suffix Optional suffix for the file name (default "").
#' @param signif.digits If not NULL, round all output values to this number of significant digits.
#'
#' @return (Invisibly) The file path to the written NetCDF file.
#'
#' @details
#' All variables and spatial grids must match the template file.
#' Attributes and dimensions are copied from the template.
#'
#' @examples
#' \dontrun{
#' # Suppose you have a template NetCDF "template.nc" with x=2, y=2, time=3
#' # Simulate some data:
#' sim_data <- replicate(4, list(
#'   precip = rnorm(3, 10, 2),
#'   temp   = rnorm(3, 20, 1)
#' ), simplify = FALSE)
#'
#' grid <- data.frame(xind = rep(1:2, each = 2), yind = rep(1:2, times = 2))
#' writeNetcdf(
#'   data = sim_data,
#'   coord.grid = grid,
#'   output.path = tempdir(),
#'   origin.date = "2000-01-01",
#'   nc.template.file = "template.nc",
#'   nc.file.prefix = "test_output"
#' )
#' }
#'
#' @importFrom ncdf4 nc_open ncvar_get ncatt_get ncvar_def ncdim_def nc_create ncvar_put ncatt_put nc_close
#' @export
writeNetcdf <- function(
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
    signif.digits = NULL) {
  # Create the main directory if not already existing
  if (!dir.exists(output.path)) {
    dir.create(output.path)
  }

  # Climate variables to write to netcdf
  variables <- names(data[[1]])

  # Simulation length
  sim.length <- nrow(data[[1]])

  # Open template netcdf file
  ncin <- ncdf4::nc_open(nc.template.file)
  on.exit(ncdf4::nc_close(ncin), add = TRUE)

  # Get dimensions from the template file
  ncin_dim <- lapply(
    1:length(names(ncin$dim)),
    function(x) ncdf4::ncvar_get(ncin, names(ncin$dim)[x])
  )
  names(ncin_dim) <- names(ncin$dim)

  # Get spatial variables and global attributes
  ncin_var_spref <- ncdf4::ncvar_get(ncin, nc.spatial.ref)
  ncin_attribs_spref <- ncdf4::ncatt_get(ncin, nc.spatial.ref)
  ncin_attribs_glob <- ncatt_get(ncin, 0)

  # Dimension names ordered as x, y, time
  ncin_dimnames <- list(
    x = ncin_attribs_spref$x_dim, y = ncin_attribs_spref$y_dim,
    time = "time"
  )

  # Set dimension variables (time, y=lat, x= long)
  ncout_dim <- list()
  ncout_dim[[ncin_dimnames$time]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$time,
    units = paste0("days since ", format(as.Date(origin.date), "%Y-%m-%d 00:00:00")),
    vals = as.numeric(1:sim.length) - 1,
    calendar = calendar.type,
    longname = NULL
  )

  ncout_dim[[ncin_dimnames$y]] <- ncdf4::ncdim_def(ncin_dimnames$y,
    units = "", vals = as.numeric(ncin_dim[[ncin_dimnames$y]]), longname = NULL
  )
  ncout_dim[[ncin_dimnames$x]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$x,
    units = "", vals = as.numeric(ncin_dim[[ncin_dimnames$x]]), longname = NULL
  )

  # Other nc attributes
  ncout_chunks <- c(
    sim.length, length(ncin_dim[[ncin_dimnames$y]]),
    length(ncin_dim[[ncin_dimnames$x]])
  )

  # Define output variables
  ncout_vars <- lapply(
    1:length(variables),
    function(x) {
      ncdf4::ncvar_def(
        name = variables[x],
        units = "",
        dim = ncout_dim,
        missval = NA,
        longname = variables[x],
        prec = "float",
        shuffle = FALSE,
        compression = nc.compression,
        chunksizes = ncout_chunks
      )
    }
  )

  # Append spatial reference variable
  ncout_vars[[length(variables) + 1]] <- ncin$var[[nc.spatial.ref]]
  ncout_vars[[length(variables) + 1]]$prec <- "integer"
  names(ncout_vars) <- c(variables, "spatial_ref")

  # template to store data from wg variables
  df0 <- array(NA, c(sim.length, ncout_dim[[ncin_dimnames$y]]$len, ncout_dim[[ncin_dimnames$x]]$len))

  # Loop through each variable and write data to netcdf
  ncout_vardata <- df0

  # Output file name
  nc_file_path <- file.path(
    output.path,
    paste0(nc.file.prefix, ifelse(nc.file.suffix != "", paste0("_", nc.file.suffix), ""), ".nc")
  )

  # create netCDF file and put arrays
  ncout_file <- ncdf4::nc_create(nc_file_path, ncout_vars, force_v4 = TRUE)
  on.exit(ncdf4::nc_close(ncout_file), add = TRUE)

  if (!is.null(signif.digits)) {
    data <- lapply(data, function(x) round(x, signif.digits))
  }

  # Loop through each variable and write data to netcdf
  for (i in 1:length(variables)) {
    ncout_vardata[1:length(ncout_vardata)] <- df0

    for (c in 1:nrow(coord.grid)) {
      ncout_vardata[, coord.grid$yind[c], coord.grid$xind[c]] <- data[[c]][[variables[i]]]
    }

    # Put variables
    ncdf4::ncvar_put(ncout_file, varid = variables[i], vals = ncout_vardata)
  }

  # Put spatial_def variable and its attributes
  ncdf4::ncvar_put(ncout_file, varid = nc.spatial.ref, vals = ncin_var_spref)
  sapply(1:length(ncin_attribs_spref), function(k) {
    ncdf4::ncatt_put(ncout_file,
      varid = nc.spatial.ref,
      attname = names(ncin_attribs_spref)[k], attval = ncin_attribs_spref[[k]]
    )
  })

  # Put global attributes
  if (length(ncin_attribs_glob) > 0) {
    sapply(
      1:length(ncin_attribs_glob),
      function(k) {
        ncdf4::ncatt_put(ncout_file,
          varid = 0,
          attname = names(ncin_attribs_glob)[k],
          attval = ncin_attribs_glob[[k]]
        )
      }
    )
  }

  # ncdf4::nc_close(ncout_file)
  invisible(nc_file_path)
}
