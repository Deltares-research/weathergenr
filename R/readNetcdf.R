
#' Read netcdf files
#'
#' Wrapper function to read files and create list object that returns a list of weather data for each grid,
#' a table of grid coordinates and a date series as well as raw data, dimensions, variables, and attributes
#' from the netcdf file.
#'
#' @param origin.date A character string to specify the date of origin
#' @param nc.path The path of the netcdf file to be read
#' @param nc.file The name of the netcdf file to be read
#' @param nc.dimnames Desired dimension names in the output data
#' @param nc.variables Desired variable names to be read from the netcdf file
#' @param has.leap.days A logical value indicating whether the gridded data includes leap days.
#' @param omit.empty.cells A lotical value indicating whether empty cells with missing values to be removed
#'
#' @return
#' @export
#' @import ncdf4
#' @import tidyr
#' @import dplyr
#' @examples
#' \dontrun{
#' output <- readNetcdf(nc.path = system.file('extdata', package = 'gridwegen'),
#'     nc.file = "ntoum.nc", nc.dimnames = list(x = "lon", y = "lat", time = "time"),
#'     nc.variables = c("precip", "temp", "temp_min", "temp_max"),
#'     origin.date = as.Date("1981-01-01"),
#'     has.leap.days = TRUE)
#'}
readNetcdf <- function(
  nc.path = NULL,
  nc.file = NULL,
  nc.dimnames = NULL,
  nc.variables = NULL,
  origin.date = NULL,
  has.leap.days = TRUE,
  omit.empty.cells = TRUE)
  {

    # Function to get dimension names
    ncGetDimNames <- function(f, v) {
        return(unlist(lapply(f$var[[v]]$dim, function(x) return(x$name))))
      }

    # Open netcdf file
    nc_in <- ncdf4::nc_open(paste(nc.path, nc.file, sep = "/"))

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
    if(isTRUE(has.leap.days)) {
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

      if(omit.empty.cells) {

        # Remove unnecessary/empty grids from tidy data table
        grd <- which(sapply(1:nrow(grid_mat), function(x)
            !is.na(mean(grid_mat$data[[x]]$temp_min))))

        grid_mat <- grid_mat[grd,]

      }


      out <- list(data = grid_mat$data,
                  coords = grid_mat %>% mutate(id = 1:n()) %>% select(-data),
                  dates = datev,
                  dimensions = nc_dims,
                  variables = nc_vars,
                  attributes = nc_attribs,
                  rawdata = nc_var_data)

  return(out)

}
