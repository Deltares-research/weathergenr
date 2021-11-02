

#' Function to read climate data from netcdf files
#'
#' @param nc.path A character string to define the path of template netcdf file.
#' @param nc.file A character string to define the name of the template netcdf file.
#' @param nc.dimnames A list object specifying the dimension names.
#' @param variables A list object specifying variable names
#' @param origin.date A character string to specify the date of origin
#' @param leap.year placeholder
#'
#' @return
#' @export
#' @importFrom tibble as_tibble
#' @importFrom ncdf4 nc_open ncatt_get ncvar_get
#' @importFrom tidyr gather
#' @importFrom dplyr filter pull mutate slice
readNetcdf <- function(
  nc.path = NULL,
  nc.file = NULL,
  nc.dimnames = NULL,
  nc.variables = NULL,
  origin.date = NULL,
  leap.year = TRUE,
  omit.nans = TRUE)

  {

    # Function to get dimension names
    ncGetDimNames <- function(f, v) {
        return(unlist(lapply(f$var[[v]]$dim, function(x) return(x$name))))
      }

    # Open netcdf file
    nc_in <- ncdf4::nc_open(paste0(nc.path, nc.file))

    # Dimensions
    nc_dim_names <- attributes(nc_in$dim)$names
    nc_dims <- lapply(1:length(nc_dim_names), function(x) ncvar_get(nc_in,nc_dim_names[x]))
    nc_dim_order <- sapply(nc.variables, function(y) sapply(nc.dimnames, function(x) which(x == ncGetDimNames(nc_in, y))))
    names(nc_dims) <- nc_dim_names

    # Variables
    nc_spatial_ref_name <- "spatial_ref"
    nc_var_names <- c(nc.variables, nc_spatial_ref_name)
    nc_vars <- lapply(nc_var_names, function(x) nc_in$var[[x]])
    names(nc_vars) <- nc_var_names

    # Variable data arrays (x=lon,y=lat,time)
    nc_var_data <- lapply(1:length(nc.variables), function(x) ncvar_get(nc_in, nc.variables[x]))
    nc_var_data <- lapply(1:length(nc_var_data), function(x) aperm(nc_var_data[[x]], nc_dim_order[,x]))
    nc_var_data[[length(nc_var_data)+1]] <- ncvar_get(nc_in, nc_spatial_ref_name)
    names(nc_var_data) <- c(nc.variables, nc_spatial_ref_name)

    # Get variable & global attributes
    nc_attribs <- lapply(nc_var_names, function(x) ncatt_get(nc_in, x))
    nc_attribs[[length(nc_attribs)+1]] <- ncatt_get(nc_in, 0)
    names(nc_attribs) <- c(nc_var_names, "global")

    # Define date vector based on calendar type (leap vs no-leap)
    if(isTRUE(leap.year)) {
      datev <- origin.date-1 + 1:length(nc_dims[[nc.dimnames$time]])
    } else {
      ymax <- round(length(nc_dims[[nc.dimnames$time]])/365)
      datev <- tibble(date = seq.Date(origin.date, origin.date + ymax*366, by = "day")) %>%
        mutate(month = as.numeric(format(date,"%m")), day = as.numeric(format(date,"%d"))) %>%
        filter(!(month == 2 & day == 29 )) %>%
        slice(1:(ymax*365)) %>% pull(date)
    }

    temp <- matrix(0, nrow = length(nc_dims$time), ncol = length(nc.variables)) %>%
      as_tibble(.name_repair = ~nc.variables)

    # translate variables to tidy format
    grid_mat <- nc_var_data[[1]][,,1] %>%
      as_tibble(.name_repair = ~as.character(1:length(nc_dims[[nc.dimnames$y]]))) %>%
      mutate(xind = 1:n(), .before = 1) %>%
      gather(yind, data, -xind) %>%
      stats::na.omit() %>%
      mutate(yind = as.numeric(yind)) %>%
      mutate(id = 1:n(), .before = xind) %>%
      mutate(x = nc_dims[[nc.dimnames$x]][xind],
             y = nc_dims[[nc.dimnames$y]][yind], .after = yind,
             data = list(temp))

      for (n in 1:nrow(grid_mat)) {
        grid_mat$data[[n]] <- as_tibble(sapply(nc_var_data[-length(nc_var_data)],
          function(x) x[grid_mat$xind[n], grid_mat$yind[n], ]))
      }

      if(omit.nans) {

        # Remove unnecessary/empty grids from tidy data table
        grd <- which(sapply(1:nrow(grid_mat), function(x)
            !is.na(mean(grid_mat$data[[x]]$temp_min))))

        grid_mat <- grid_mat[grd,]

      }


      out <- list(tidy_data = grid_mat %>% mutate(id = 1:n()),
                  nc_dates = datev,
                  nc_dimensions = nc_dims,
                  nc_variables = nc_vars,
                  nc_variable_data = nc_var_data,
                  nc_attributes = nc_attribs)

  return(out)

}
