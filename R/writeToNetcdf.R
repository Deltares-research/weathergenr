
#' A function to write a netcdf file from tidydata
#'
#' @param in.path A character string to define the path of template netcdf file.
#' @param in.file A character string to define the name of the template netcdf file.
#' @param out.path A character string to define the path of resulting netcdf file.
#' @param out.file A character string to define the name of resulting netcdf file.
#' @param origin.date A character string to set the origin date in the resulting netcdf file.
#' @param gridded.data A multidimensional array to be written to the netcdf file.
#' @param out.variables A character vector to define which variables to be outputted to the netcdf file.
#' @param dim.names A character vector of dimension names.
#' @param verbose A logical to print more information on the screen.
#'
#' @return
#' @export
#' @import ncdf4
writeToNetcdf <- function(
  in.path = NULL,
  in.file = NULL,
  out.path = NULL,
  out.file = NULL,
  out.variables = NULL,
  dim.names = NULL,
  origin.date = NULL,
  gridded.data = NULL,
  verbose = FALSE)

  {

  # Open template netcdf file
  ncin <- nc_open(paste0(in.path, in.file))

  # Global attributes
  nc_atts <- ncatt_get(ncin, 0)


  #::::::::::::::::::::::: DIMENSIONS ::::::::::::::::::::::::::::::::::::::::::

  # Define time dimension
  time_vals  <- 1:length(wg_sim_dates)
  time_units <- paste0("days since ",
                       format(round(as.POSIXct(origin.date), units = "day"), '%Y-%m-%d %M:%H:%S'))

  dim_names <- attributes(ncin$dim)$names
  DIMS <- vector(mode = "list", length = length(dim_names))
  names(DIMS) <- dim_names

  # Set spatial dimensions
  DIMS[[dim.names$x]] <- ncdim_def(dim.names$x,
        units= "", vals = ncvar_get(ncin,dim.names$x))
  DIMS[[dim.names$y]] <- ncdim_def(dim.names$x,
        units= "", vals = ncvar_get(ncin,dim.names$y))

  # Set time dimension
  DIMS[[dim.names$time]] <- ncdim_def(dim.names$time,
      units = time_units, vals = time_vals, calendar = "noleap")


  #::::::::::::::::::::: VARIABLES :::::::::::::::::::::::::::::::::::::::::::::

  # GET ALL THE VARIABLES
  var_names <- attributes(ncin$var)$names
  VARS <- vector(mode = "list", length = length(var_names))
  ARRAY <- vector(mode = "list", length = length(var_names))
  ATTRIBS <- vector(mode = "list", length = length(var_names))
  names(VARS) <- var_names
  names(ARRAY) <- var_names
  names(ATTRIBS) <- var_names

  #Replace precision parameter for all variables
  var_prec <- lapply(var_names, function(x) ncin$var[[x]]$prec)
  var_prec <- lapply(var_prec, function(x) str_replace(x, "int", "integer"))
  var_prec <- lapply(var_prec, function(x) str_replace_all(x, "unsigned byte", "byte"))
  names(var_prec) <- var_names

  # Loop through variables
  for (i in 1:length(var_names)) {

    # Current variable
    varc <- var_names[i]

    # Set missing values based on precision
    if(var_prec[[i]] %in% c("float", "double")) {
      varc_missingValue <- NA
    } else if (var_prec[[i]] == "integer") {
      varc_missingValue <- -1
    } else {
      varc_missingValue <- NULL
    }

    # Define variable dimensions
    if(is.null(ncin$var[[var_names[i]]]$dim)) {
      varc_dims <- list()
    } else {
      varc_dims <- DIMS[sapply(ncin$var[[varc]]$dim, function(x) x$name)]
    }

    VARS[[varc]] <- ncvar_def(
      name = varc,
      units = ncin$var[[varc]]$units,
      dim = varc_dims,
      missval = varc_missingValue,
      longname = ncin$var[[varc]]$longname,
      prec =  var_prec[[varc]],
      shuffle = FALSE,
      compression = ncin$var[[varc]]$compression)

  }

  # TRANSLATE INTO TIDY-FORMAT
  coordGrid <- ncvar_get(ncin,out.variables[1])[,,1] %>% as_tibble() %>%
    setNames(1:DIMS[[dim.names$y]]$len) %>%
    mutate(xind = 1:n(), .before = 1) %>%
    gather(yind, data, -xind) %>%
    mutate(yind = as.numeric(yind)) %>%
    mutate(x = DIMS[[dim.names$x]]$vals[xind],
           y = DIMS[[dim.names$y]]$vals[yind], .after = yind) %>%
    na.omit() %>%
    mutate(data = list(NA))


  # template to store data from wg variables
  array_temp <- array(NaN, c(DIMS[[dim.names$x]]$len,
                          DIMS[[dim.names$y]]$len,
                          DIMS[[dim.names$time]]$len))

  #Loop through each variable
  for (i in 1:length(var_names)) {

    # Current variable
    varc <- var_names[i]
    if(varc %in% out.variables) {

      ARRAY[[varc]] <- array_temp

      for (x in 1:nrow(coordGrid)) {

        xind_i <- coordGrid$xind[x]
        yind_i <- coordGrid$yind[x]

        ARRAY[[varc]][xind_i,yind_i, ] <- gridded.data[[i]][[varc]]
      }

    } else {
      ARRAY[[varc]] <- ncvar_get(ncin,varc)
    }


    # Get attributes for each variable
    ATTRIBS[[varc]] <- ncatt_get(ncin,varc)
    if(is.null(ATTRIBS[[i]]$units)) {ATTRIBS[[i]]$units = "not specified"}


  }

  #Remove missing value attribute
  ATTRIBS <- lapply(1:length(VARS),
                    function(x) within(ATTRIBS[[x]], rm(`_FillValue`)))


  # ::::::::: CREATE NETCDF AND PUT VARIABLES ::::::::::::::::::::::::::::::::::

  # create netCDF file and put arrays
  ncout <- nc_create(paste0(out.path, out.file), VARS, force_v4 = TRUE)
  for (i in 1:length(var_names)) {

    # Current variable
    varc <- var_names[i]

    # Put variables
    ncvar_put(ncout, varid = VARS[[varc]], vals = ARRAY[[varc]])

    # Put additional attributes for each variable
    if(!is.null(ATTRIBS[[varc]])) {

      sapply(1:length(ATTRIBS[[varc]]), function(k)
        ncatt_put(ncout, varid = varc, attname = names(ATTRIBS[[varc]])[k],
                  attval = ATTRIBS[[varc]][[k]]))
    }
  }
  # put global attributes
  if(length(nc_atts)>0) {
    sapply(1:length(nc_atts), function(k)
      ncatt_put(ncout, varid = 0, attname = names(nc_atts)[k],
                attval = nc_atts[[k]]))
  }

  nc_close(ncout)
  if (isTRUE(verbose)) {
    print(paste0(out.file, " created successfully!"))
  }


}



