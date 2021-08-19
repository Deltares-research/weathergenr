

#' Function to read climate data from netcdf files
#'
#' @param in.path A character string to define the path of template netcdf file.
#' @param in.file A character string to define the name of the template netcdf file.
#' @param dim.names A list object specifying the dimension names.
#' @param variables A list object specifying variable names
#' @param origin.date A character string to specify the date of origin
#'
#' @return
#' @export
#' @import ncdf4
#' @import tidyr
#' @import dplyr
readNetcdf <- function(
  in.path = NULL,
  in.file = NULL,
  dim.names = NULL,
  variables = NULL,
  origin.date = NULL,
  leap.year = TRUE)
  {

   ncGetDimNames <- function(f, v) {
      return(unlist(lapply(f$var[[v]]$dim, function(x) return(x$name))))
    }

  # ::::::::: Load & prepare climate data ::::::::::::::::::::::::::::::::::::::

  # Open template netcdf file
  nc_in <- ncdf4::nc_open(paste0(in.path, in.file))

  # GET DIMENSIONS
  nc_dim_names <- attributes(nc_in$dim)$names
  nc_dims <- lapply(1:length(nc_dim_names), function(x) ncvar_get(nc_in,nc_dim_names[x]))
  nc_dim_order <- sapply(variables, function(y) sapply(dim.names, function(x) which(x == ncGetDimNames(nc_in, y))))
  names(nc_dims) <- nc_dim_names

  # GET VARIABLES
  nc_spatial_ref_name <- "spatial_ref"
  nc_var_names <- c(variables, nc_spatial_ref_name)
  nc_vars <- lapply(nc_var_names, function(x) nc_in$var[[x]])
  names(nc_vars) <- nc_var_names

  # GET VARIABLE DATA ARRAYS AND REORDER DIMENSIONS (x=lon,y=lat,time)
  nc_var_data <- lapply(1:length(variables), function(x) ncvar_get(nc_in, variables[x]))
  nc_var_data <- lapply(1:length(nc_var_data), function(x) aperm(nc_var_data[[x]], nc_dim_order[,x]))
  nc_var_data[[length(nc_var_data)+1]] <- ncvar_get(nc_in, nc_spatial_ref_name)
  names(nc_var_data) <- c(variables, nc_spatial_ref_name)

  # GET ATTRIBUTES
  nc_attribs <- lapply(nc_var_names, function(x) ncatt_get(nc_in, x))
  nc_attribs[[length(nc_attribs)+1]] <- ncatt_get(nc_in, 0)
  names(nc_attribs) <- c(nc_var_names, "global")

  # Define date vector based on leap_day choice
  if(isTRUE(leap.year)) {
    datev <- origin.date-1 + 1:length(nc_dims[[dim.names$time]])
  } else {
    ymax <- round(length(nc_dims[[dim.names$time]])/365)
    datev <- tibble(date = seq.Date(origin.date, origin.date + years(ymax), by = "day")) %>%
      filter(!(month(date) == 2 & (day(date) == 29))) %>% slice(1:(ymax*365)) %>% pull(date)
  }


  temp <- matrix(0, nrow = length(nc_dims$time), ncol = length(variables)) %>%
    as_tibble(.name_repair = ~variables) %>%
    mutate(date = datev, .before = 1)

  # translate variables to tidy format
  grid_mat <- nc_var_data[[1]][,,1] %>%
    as_tibble(.name_repair = ~as.character(1:length(nc_dims[[dim.names$y]]))) %>%
    mutate(xind = 1:n(), .before = 1) %>%
    gather(yind, data, -xind) %>%
    na.omit() %>%
    mutate(yind = as.numeric(yind)) %>%
    mutate(id = 1:n(), .before = xind) %>%
    mutate(x = nc_dims[[dim.names$x]][xind],
           y = nc_dims[[dim.names$y]][yind], .after = yind,
           data = list(temp))

  for (n in 1:nrow(grid_mat)) {
    grid_mat$data[[n]][,-1] <- sapply(nc_var_data[-length(nc_var_data)],
      function(x) x[grid_mat$xind[n], grid_mat$yind[n], ])
  }

    out <- list(tidy_data = grid_mat,
                nc_dimensions = nc_dims,
                nc_variables = nc_vars,
                nc_variable_data = nc_var_data,
                nc_attributes = nc_attribs)



  # GET SPATIAL REF VARIABLE AND ATTRIBUTES

  # nc_spatial_ref_value <- ncvar_get(nc_in, nc_spatial_ref_name)
  # nc_spatial_ref_attribs <- ncatt_get(nc_in,nc_spatial_ref_name)
  # nc_spatial_ref_def <- ncvar_def(name = nc_spatial_ref_name,
  #   units = nc_in$var[[nc_spatial_ref_name]]$units,
  #   dim = nc_in$var[[nc_spatial_ref_name]]$dims, prec = "integer")
  #


  return(out)

}
