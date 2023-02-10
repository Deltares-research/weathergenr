

#' Write to netcdf file
#'
#' @param output.path Placeholder
#' @param origin.date Placeholder
#' @param calendar.type Placeholder
#' @param nc.compression Placeholder
#' @param nc.template.file placeholder
#' @param data placeholder
#' @param coord.grid placeholder
#' @param nc.file.prefix placeholder
#' @param nc.file.suffix placeholder
#' @param nc.spatial.ref placeholder
#'
#' @return
#' @export
#'
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
  nc.file.suffix = "")

{

  # Create the main directory if not already existing
  if (!dir.exists(output.path)) {dir.create(output.path)}

  # Climate variables to write to netcdf
  variables <- names(data[[1]])

  # Simulation length
  sim.length = nrow(data[[1]])

  # Open netcdf file
  ncin <- ncdf4::nc_open(nc.template.file)

  # Get dimension variables
  ncin_dim <- lapply(1:length(names(ncin$dim)), function(x) ncdf4::ncvar_get(ncin,names(ncin$dim)[x]))
  names(ncin_dim) <- names(ncin$dim)

  # Get variables & attributes
  ncin_var_spref <- ncdf4::ncvar_get(ncin, nc.spatial.ref)
  ncin_attribs_spref <-ncdf4:: ncatt_get(ncin, nc.spatial.ref)
  ncin_attribs_glob <- ncatt_get(ncin, 0)

  # Dimension names ordered as x, y, time
  ncin_dimnames <- list(x = ncin_attribs_spref$x_dim, y = ncin_attribs_spref$y_dim,
    time = "time")

  # Set dimension variables (time, y=lat, x= long)
  ncout_dim <- list()
  ncout_dim[[ncin_dimnames$time]] <- ncdf4::ncdim_def(name = ncin_dimnames$time,
      units = paste0("days since ", format(round(as.POSIXct(origin.date)), '%Y-%m-%d %M:%H:%S')),
      vals = 1:sim.length, calendar = calendar.type)
  ncout_dim[[ncin_dimnames$y]] <- ncdf4::ncdim_def(ncin_dimnames$y,
        units= "", vals = ncin_dim[[ncin_dimnames$y]])
  ncout_dim[[ncin_dimnames$x]] <- ncdf4::ncdim_def(name = ncin_dimnames$x,
        units= "", vals = ncin_dim[[ncin_dimnames$x]])

  # Other nc attributes
  ncout_chunks <- c(1, length(ncin_dim[[ncin_dimnames$y]]), length(ncin_dim[[ncin_dimnames$x]]))

  # Define output variables
  ncout_vars <- lapply(1:length(variables),
    function(x) ncdf4::ncvar_def(
        name = variables[x],
        units = "",
        dim = ncout_dim,
        missval = NA,
        longname = variables[x],
        prec =  "float",
        shuffle = FALSE,
        compression = nc.compression,
        chunksizes= ncout_chunks)
  )

  # Append spatial reference variable
  ncout_vars[[length(ncout_vars)+1]] <- ncin$var[[nc.spatial.ref]]
  #ncout_vars[length(ncout_vars)][[1]]$prec <- "integer"
  names(ncout_vars) <- c(variables, "spatial_ref")

  # template to store data from wg variables
  df0 <- array(NA, c(ncout_dim[[ncin_dimnames$time]]$len,
    ncout_dim[[ncin_dimnames$y]]$len, ncout_dim[[ncin_dimnames$x]]$len))

  #Loop through each variable and write data to netcdf
  ncout_vardata <- df0

  # create netCDF file and put arrays
  ncout_file <- ncdf4::nc_create(file.path(output.path, nc.file.prefix,"_",nc.file.suffix, ".nc"),
    ncout_vars, force_v4 = TRUE)

  #Loop through each variable and write data to netcdf
  for (i in 1:length(variables)) {

    ncout_vardata[1:length(ncout_vardata)] <- df0

    for (c in 1:nrow(coord.grid)) {
      ncout_vardata[, coord.grid$yind[c], coord.grid$xind[c]] <- data[[c]][[variables[i]]]
    }

    # Put variables
    ncdf4::ncvar_put(ncout_file, varid = variables[i], vals = ncout_vardata)

  }

  # Put spatial_def variable and attributes
  ncdf4::ncvar_put(ncout_file, varid = nc.spatial.ref, vals = ncin_var_spref)
  sapply(1:length(ncin_attribs_spref), function(k)
    ncdf4::ncatt_put(ncout_file, varid = nc.spatial.ref,
              attname = names(ncin_attribs_spref)[k], attval = ncin_attribs_spref[[k]]))

 # Put global attributes
  if(length(ncin_attribs_glob)>0) {

    sapply(1:length(ncin_attribs_glob),
      function(k) ncdf4::ncatt_put(ncout_file, varid = 0,
              attname = names(ncin_attribs_glob)[k],
              attval = ncin_attribs_glob[[k]]))
  }

  ncdf4::nc_close(ncout_file)

}


