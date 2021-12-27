
#' Read netcdf files
#'
#' Wrapper function to read files and create list object that returns a list of weather data for each grid,
#' a table of grid coordinates and a date series as well as raw data, dimensions, variables, and attributes
#' from the netcdf file.
#'
#' @param nc.file The name of the netcdf file to be read
#' @param has.leap.days A logical value indicating whether the gridded data includes leap days.
#' @param omit.empty.cells A lotical value indicating whether empty cells with missing values to be removed
#'
#' @return
#' @export
#' @import ncdf4
#' @import tidyr
#' @import dplyr
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#' output <- readNetcdf(nc.path = system.file('extdata', package = 'gridwegen'),
#'     nc.file = "ntoum.nc", nc_dimnames = list(x = "lon", y = "lat", time = "time"),
#'     nc.variables = c("precip", "temp", "temp_min", "temp_max"),
#'     origin_date = as.Date("1981-01-01"),
#'     has.leap.days = TRUE)
#'}
readNetcdf <- function(
  nc.file = NULL,
  has.leap.days = TRUE,
  omit.empty.cells = TRUE,
  spatial.ref.variable = "spatial_ref")

  {

    #Workaround for rlang warning
    month <- day <- xind <- yind <- month <- day <- wyear <- 0

    # function to get dimension names
    ncGetDimNames <- function(f, v) {
        return(unlist(lapply(f$var[[v]]$dim, function(x) return(x$name))))
    }

    # Open netcdf file
    nc_in <- ncdf4::nc_open(nc.file)

    # Get dimension variables
    nc_dim_names <- attributes(nc_in$dim)$names
    nc_dim <- lapply(1:length(nc_dim_names), function(x) ncvar_get(nc_in,nc_dim_names[x]))
    names(nc_dim) <- nc_dim_names

    # Get spatial reference variable
    nc_var_data_spref <- ncvar_get(nc_in, spatial.ref.variable)
    nc_attribs_spref <- ncatt_get(nc_in, spatial.ref.variable)

    time_unit <- ncatt_get(nc_in, "time")$units
    origin_date <- as.Date(sub("^\\D+", "", time_unit))


    # Dimension names ordered as x, y, time
    nc_dimnames <- list(x = nc_attribs_spref$x_dim,
       y = nc_attribs_spref$y_dim,time = "time")

    # Get variables
    nc_var_names <- setdiff(attributes(nc_in$var)$names, spatial.ref.variable)
    nc_var_data <- lapply(nc_var_names, function(x) ncvar_get(nc_in, x))
    nc_dimord <- sapply(nc_var_names, function(y)
      sapply(nc_dimnames, function(x) which(x == ncGetDimNames(nc_in, y))))
    nc_var_data <- lapply(1:length(nc_var_data), function(x) aperm(nc_var_data[[x]], nc_dimord[,x]))
    names(nc_var_data) <- nc_var_names

    # Get local and global attributes
    nc_attribs <- lapply(names(nc_var_data), function(x) ncatt_get(nc_in, x))
    names(nc_attribs) <- nc_var_names
    nc_attribs_glob <- ncatt_get(nc_in, 0)

    # Define date vector based on calendar type
    if(isTRUE(has.leap.days)) {
      datev <- origin_date-1 + 1:length(nc_dim[[nc_dimnames$time]])
    } else {
      ymax <- round(length(nc_dim[[nc_dimnames$time]])/365)
      datev <- tibble(date = seq.Date(origin_date, origin_date + ymax*366, by = "day")) %>%
        mutate(month = as.numeric(format(date,"%m")), day = as.numeric(format(date,"%d"))) %>%
        filter(!(month == 2 & day == 29 )) %>%
        slice(1:(ymax*365)) %>% pull(date)
    }

    # Template data table for grid
    df0 <- matrix(0, nrow = length(nc_dim$time), ncol = length(nc_var_names)) %>%
      as_tibble(.name_repair = ~nc_var_names)

    # translate variables to tidy format
    grid_mat <- nc_var_data[[1]][,,1] %>%
      as_tibble(.name_repair = ~as.character(1:length(nc_dim[[nc_dimnames$y]]))) %>%
      mutate(xind = 1:n(), .before = 1) %>%
      gather(yind, data, -xind) %>%
      stats::na.omit() %>%
      mutate(yind = as.numeric(yind)) %>%
      mutate(id = 1:n(), .before = xind) %>%
      mutate(x = nc_dim[[nc_dimnames$x]][xind],
             y = nc_dim[[nc_dimnames$y]][yind], .after = yind, data = list(df0))


      for (n in 1:nrow(grid_mat)) {
        grid_mat$data[[n]] <- as_tibble(sapply(nc_var_data,
          function(x) x[grid_mat$xind[n], grid_mat$yind[n], ]))
      }

      if(omit.empty.cells) {

        # Remove unnecessary/empty grids from tidy data table
        grd <- which(sapply(1:nrow(grid_mat), function(x) !is.na(sum(grid_mat$data[[x]]))))
        grid_mat <- grid_mat[grd,]

      }

      out <- list(data = grid_mat$data,
                  coords = grid_mat %>% mutate(id = 1:n()) %>% select(-data),
                  dates = datev,
                  dimensions = nc_dim,
                  attributes = c(nc_attribs, setNames(list(nc_attribs_spref),
                    spatial.ref.variable), global = list(nc_attribs_glob)),
                  rawdata = nc_var_data)

  return(out)

}
