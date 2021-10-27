

#' Title
#'
#' @param output.path Placeholder
#' @param nc.dimensions Placeholder
#' @param nc.dimnames Placeholder
#' @param origin.date Placeholder
#' @param sim.length Placeholder
#' @param calendar.type Placeholder
#' @param variables Placeholder
#' @param variable.units Placeholder
#' @param nc.compression Placeholder
#'
#' @return
#' @export
#'
writeNetcdf <- function(
  data = NULL,
  output.path = NULL,
  nc.dimensions = NULL,
  nc.dimnames = NULL,
  origin.date = NULL,
  calendar.type = "no leap",
  variables = NULL,
  variable.units = NULL,
  nc.compression = 4,
  file.suffix = 1

) {

  if (!dir.exists(output.path)) {dir.create(output.path)}

  sim.length = nrow(data[[1]])

  # Dimensions in the outputted variable (order matters!)
  dim_ord <- nc.dimensions

  ncout_dims <- list()
  ncout_dims[[nc.dimnames$time]] <- ncdf4::ncdim_def(name = nc.dimnames$time,
      units = paste0("days since ", format(round(as.POSIXct(origin.date), units = "day"),
        '%Y-%m-%d %M:%H:%S')), vals = 1:sim.length, calendar = "no leap")

  ncout_dims[[nc.dimnames$y]] <- ncdf4::ncdim_def(nc.dimnames$y,
        units= "", vals = nc_data$nc_dimensions[[nc.dimnames$y]])

  ncout_dims[[nc.dimnames$x]] <- ncdf4::ncdim_def(name = nc.dimnames$x,
        units= "", vals = nc_data$nc_dimensions[[nc.dimnames$x]])

  ncout_dim_reorder <- unname(sapply(names(nc.dimnames), function(x) which(x == dim_ord)))

  # Other nc attributes
  ncout_chunksize <- c(1,length(nc_data$nc_dimensions[[2]]),length(nc_data$nc_dimensions[[3]]))

  # ncout variables
  ncout_vars <- lapply(1:length(variables), function(x) ncdf4::ncvar_def(
        name = variables[x],
        units = variable.units[x],
        dim = ncout_dims,
        missval = NA,
        longname = variables[x],
        prec =  "float",
        shuffle = FALSE,
        compression = nc.compression,
        chunksizes= ncout_chunksize))

  ncout_vars[[length(ncout_vars)+1]] <- nc_data$nc_variables$spatial_ref
  names(ncout_vars) <- c(variables, "spatial_ref")

  # TRANSLATE INTO TIDY-FORMAT
  coordGrid <- climate_tidy %>% mutate(date = list(NA))

  # template to store data from wg variables
  var_empty <- array(NA, c(
    ncout_dims[[nc.dimnames$time]]$len,
    ncout_dims[[nc.dimnames$y]]$len,
    ncout_dims[[nc.dimnames$x]]$len))

  #Loop through each variable and write data to netcdf
  ncout_vardata <- var_empty

  # create netCDF file and put arrays
  ncout_file <- ncdf4::nc_create(paste0(output.path, "hist_rlz_",file.suffix, ".nc"),
    ncout_vars, force_v4 = TRUE)

  #Loop through each variable and write data to netcdf
  for (i in 1:length(variables)) {

    ncout_vardata[1:length(ncout_vardata)] <- var_empty

    for (c in 1:nrow(coordGrid)) {
      ncout_vardata[, coordGrid$yind[c], coordGrid$xind[c]] <- data[[c]][[variables[i]]]
    }

    # Put variables
    ncdf4::ncvar_put(ncout_file, varid = variables[i], vals = ncout_vardata)

  }

  # Put spatial_def variable and attributes
  ncdf4::ncvar_put(ncout_file, varid = nc_data$nc_variables$spatial_ref$name,
            vals = nc_data$nc_variable_data$spatial_ref)

  sapply(1:length(nc_data$nc_attributes$spatial_ref), function(k)
    ncdf4::ncatt_put(ncout_file,
              varid = nc_data$nc_variables$spatial_ref$name,
              attname = names(nc_data$nc_attributes$spatial_ref)[k],
              attval = nc_data$nc_attributes$spatial_ref[[k]]))

 # Put global attributes
  if(length(nc_data$nc_attributes$global)>0) {

    sapply(1:length(nc_data$nc_attributes$global),
      function(k) ncdf4::ncatt_put(ncout_file, varid = 0,
              attname = names(nc_data$nc_attributes$global)[k],
              attval = nc_data$nc_attributes$global[[k]]))
  }

  ncdf4::nc_close(ncout_file)

}


