
#' Read netcdf files
#'
#' Wrapper function to read files and create list object that returns a list of weather data for each grid,
#' a table of grid coordinates and a date series as well as raw data, dimensions, variables, and attributes
#' from the netcdf file.
#'
#' @param nc.file The name of the netcdf file to be read
#' @param leap.days A logical value indicating whether the gridded data includes leap days.
#' @param omit.empty A lotical value indicating whether empty cells with missing values to be removed
#' @param spatial.ref PLACEHOLDER
#' @param signif.digits PLACEHOLDER
#'
#' @return
#' @export
#' @import ncdf4
#' @import tidyr
#' @import dplyr
#' @examples
#' \dontrun{
#' output <- readNetcdf(nc.path = system.file('extdata', package = 'gridwegen'),
#'     nc.file = "ntoum.nc", nc_dimnames = list(x = "lon", y = "lat", time = "time"),
#'     nc.variables = c("precip", "temp", "temp_min", "temp_max"),
#'     origin_date = as.Date("1981-01-01"),
#'     leap.days = TRUE)
#'}
readNetcdf <- function(
  nc.file = NULL,
  leap.days = TRUE,
  omit.empty = TRUE,
  spatial.ref = "spatial_ref",
  signif.digits = NULL)

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
    nc_dim <- lapply(1:length(nc_dim_names), function(x) ncdf4::ncvar_get(nc_in,nc_dim_names[x]))
    names(nc_dim) <- nc_dim_names

    # Get spatial reference variable
    nc_var_data_spref <- ncvar_get(nc_in, spatial.ref)
    nc_attribs_spref <- ncatt_get(nc_in, spatial.ref)

    # Get time parameters
    time_unit <- ncatt_get(nc_in, "time")$units
    origin_date <- as.Date(sub("^\\D+", "", time_unit))

    # Dimension names ordered as x, y, time
    nc_dimnames <- list(x = nc_attribs_spref$x_dim, y = nc_attribs_spref$y_dim, time = "time")

    # Get variables
    nc_var_names <- setdiff(attributes(nc_in$var)$names, spatial.ref)
    nc_var_data <- sapply(nc_var_names, function(x) ncdf4::ncvar_get(nc_in, x), simplify = "array")
    if(!is.null(signif.digits)) nc_var_data <- round(nc_var_data, signif.digits)

    nc_dimord <- sapply(nc_var_names, function(y)
      sapply(nc_dimnames, function(x) which(x == ncGetDimNames(nc_in, y))))
    nc_dimord <- rbind(nc_dimord , rep(4, length(nc_var_names)))
    nc_var_data <- aperm(nc_var_data, nc_dimord[,1])
    names(nc_var_data) <- nc_var_names

    # Get local and global attributes
    nc_attribs <- lapply(nc_var_names, function(x) ncatt_get(nc_in, x))
    names(nc_attribs) <- nc_var_names
    nc_attribs_glob <- ncatt_get(nc_in, 0)

    #nc_var_names <- setdiff(attributes(nc_in$var)$names, spatial.ref)
    #nc_var_data <- lapply(nc_var_names, function(x) ncdf4::ncvar_get(nc_in, x))
    #nc_dimord <- sapply(nc_var_names, function(y)
    #  sapply(nc_dimnames, function(x) which(x == ncGetDimNames(nc_in, y))))
    #nc_var_data <- lapply(1:length(nc_var_data), function(x) aperm(nc_var_data[[x]], nc_dimord[,x]))
    #names(nc_var_data) <- nc_var_names

    # Get local and global attributes
    #nc_attribs <- lapply(names(nc_var_data), function(x) ncatt_get(nc_in, x))
    #names(nc_attribs) <- nc_var_names
    #nc_attribs_glob <- ncatt_get(nc_in, 0)

    # Define date vector based on calendar type
    if(isTRUE(leap.days)) {
      datev <- origin_date + nc_dim[[nc_dimnames$time]]
    } else {
      ymax <- round(length(nc_dim[[nc_dimnames$time]])/365)
      begin_date <- origin_date + nc_dim[[nc_dimnames$time]][[1]]
      datev <- tibble(date = seq.Date(begin_date, begin_date + ymax*366, by = "day")) %>%
        mutate(month = as.numeric(format(date,"%m")), day = as.numeric(format(date,"%d"))) %>%
        filter(!(month == 2 & day == 29 )) %>%
        slice(1:(ymax*365)) %>% pull(date)
    }

    # translate variables to tidy format
    grid_mat <- nc_var_data[,,1,1] %>%
      as_tibble(.name_repair = ~as.character(1:length(nc_dim[[nc_dimnames$y]]))) %>%
      mutate(xind = 1:n(), .before = 1) %>%
      gather(yind, data, -xind) %>%
      stats::na.omit() %>%
      mutate(yind = as.numeric(yind)) %>%
      mutate(id = 1:n(), .before = xind) %>%
      mutate(x = nc_dim[[nc_dimnames$x]][xind],
             y = nc_dim[[nc_dimnames$y]][yind], .after = yind) %>%
      select(-data)

    #grid_mat <- nc_var_data[[1]][,,1] %>%
    #  as_tibble(.name_repair = ~as.character(1:length(nc_dim[[nc_dimnames$y]]))) %>%
    #  mutate(xind = 1:n(), .before = 1) %>%
    #  gather(yind, data, -xind) %>%
    #  stats::na.omit() %>%
    #  mutate(yind = as.numeric(yind)) %>%
    #  mutate(id = 1:n(), .before = xind) %>%
    #  mutate(x = nc_dim[[nc_dimnames$x]][xind],
    #         y = nc_dim[[nc_dimnames$y]][yind], .after = yind)

    # Write gridded data to tidy data frames
    df0 <- matrix(0, nrow = length(nc_dim$time), ncol = length(nc_var_names)) %>%
      as_tibble(.name_repair = ~nc_var_names)
    data <- rep(list(df0), nrow(grid_mat))

    x_iter <- grid_mat$xind
    y_iter <- grid_mat$yind

    for (n in 1:nrow(grid_mat)) {
      data[[n]][] <- nc_var_data[x_iter[n], y_iter[n], ,]
    }

    # Write gridded data to tidy data frames
    #df0 <- matrix(0, nrow = length(nc_dim$time), ncol = length(nc_var_names)) %>%
    #  as_tibble(.name_repair = ~nc_var_names)
    #data <- rep(list(df0), nrow(grid_mat))
    #for (n in 1:nrow(grid_mat)) {
    #  data[[n]] <- as_tibble(sapply(nc_var_data,
    #      function(x) x[grid_mat$xind[n], grid_mat$yind[n], ]))
    #}

    # Remove unnecessary/empty grids from tidy data table
    if(omit.empty) {
      grd <- which(sapply(1:nrow(grid_mat), function(x) !is.na(sum(data[[x]]))))
      grid_mat <- grid_mat[grd,]
    }

    return(list(data = data,
                grid = grid_mat %>% mutate(id = 1:n()),
                date = datev,
                dimensions = nc_dim,
                attributes = c(nc_attribs, stats::setNames(list(nc_attribs_spref),
                               spatial.ref), global = list(nc_attribs_glob))))
}
