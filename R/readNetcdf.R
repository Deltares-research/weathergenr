#' Read a gridded NetCDF file and return weather data in tidy format
#'
#' @param nc.file Path to the NetCDF file.
#' @param leap.days Logical: does data include leap days? (default: TRUE)
#' @param omit.empty Logical: remove empty/missing-value cells? (default: TRUE)
#' @param spatial.ref String: name of spatial reference variable in NetCDF (default: "spatial_ref")
#' @param signif.digits (optional) Number of significant digits to round data.
#'
#' @return List with elements:
#'   \item{data}{List of tibbles of weather variables for each grid cell (length = number of cells)}
#'   \item{grid}{Tibble of grid coordinates and indices}
#'   \item{date}{Vector of Date objects}
#'   \item{dimensions}{Named list of dimension variables (time, x, y, etc.)}
#'   \item{attributes}{List of attributes for variables and global attributes}
#'
#' @import ncdf4
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @export
readNetcdf <- function(
  nc.file,
  leap.days = TRUE,
  omit.empty = TRUE,
  spatial.ref = "spatial_ref",
  signif.digits = NULL
) {

    # Validate file
    if (!file.exists(nc.file)) stop("File does not exist: ", nc.file)

    #Workaround for rlang warning
    month <- day <- xind <- yind <- month <- day <- wyear <- 0

    # function to get dimension names
    ncGetDimNames <- function(f, v) {
        return(unlist(lapply(f$var[[v]]$dim, function(x) return(x$name))))
    }

    # Open netcdf file
    nc_in <- ncdf4::nc_open(nc.file)
    on.exit(ncdf4::nc_close(nc_in))

    # Gather dimension names and variables
    nc_dim_names <- attributes(nc_in$dim)$names
    nc_dim <- lapply(1:length(nc_dim_names), function(x) ncdf4::ncvar_get(nc_in,nc_dim_names[x]))
    names(nc_dim) <- nc_dim_names

    # Read spatial reference (if present)
    if (spatial.ref %in% names(nc_in$var)) {
      nc_attribs_spref <- ncatt_get(nc_in, spatial.ref)
      nc_var_data_spref <- ncvar_get(nc_in, spatial.ref)
    } else {
      stop("Spatial reference variable '", spatial.ref, "' not found.")
    }

    # Read time units and origin date
    time_unit <- ncatt_get(nc_in, "time")$units
    origin_date <- tryCatch(
      as.Date(sub("^\\D+", "", time_unit)),
      error = function(e) stop("Could not extract origin date from NetCDF time units.")
    )

    # Dimension names ordered as x, y, time
    nc_dimnames <- list(x = nc_attribs_spref$x_dim, y = nc_attribs_spref$y_dim, time = "time")

    # Get variable names (excluding spatial reference if present)
    nc_var_names <- setdiff(attributes(nc_in$var)$names, spatial.ref)
    nc_var_data <- sapply(nc_var_names, function(x) ncdf4::ncvar_get(nc_in, x), simplify = "array")
    if(!is.null(signif.digits)) nc_var_data <- round(nc_var_data, signif.digits)

    nc_dimord <- sapply(nc_var_names, function(y)
      sapply(nc_dimnames, function(x) which(x == ncGetDimNames(nc_in, y))))
    nc_dimord <- rbind(nc_dimord , rep(4, length(nc_var_names)))
    nc_var_data <- aperm(nc_var_data, nc_dimord[,1])
    names(nc_var_data) <- nc_var_names

    # Get variable and global attributes
    nc_attribs <- lapply(nc_var_names, function(x) ncatt_get(nc_in, x))
    names(nc_attribs) <- nc_var_names
    nc_attribs_glob <- ncatt_get(nc_in, 0)

    # Define date vector based on calendar type
    time_dim <- nc_dim[[nc_dimnames$time]]
    if (is.null(time_dim)) stop("'time' dimension not found in NetCDF.")
    if(leap.days) {
      datev <- origin_date + time_dim
    } else {
      # Remove leap days
      ymax <- round(length(time_dim)/365)
      begin_date <- origin_date + time_dim[[1]]
      datev <- tibble::tibble(date = seq.Date(begin_date, begin_date + ymax*366, by = "day")) %>%
        dplyr::mutate(month = as.numeric(format(date,"%m")),
                      day = as.numeric(format(date,"%d"))) %>%
        dplyr::filter(!(month == 2 & day == 29 )) %>%
        dplyr::slice(1:(ymax*365)) %>%
        dplyr::pull(date)
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

    # Write gridded data to tidy data frames
    df0 <- matrix(0, nrow = length(nc_dim$time), ncol = length(nc_var_names)) %>%
      as_tibble(.name_repair = ~nc_var_names)
    data <- rep(list(df0), nrow(grid_mat))

    x_iter <- grid_mat$xind
    y_iter <- grid_mat$yind

    for (n in 1:nrow(grid_mat)) {
      data[[n]][] <- nc_var_data[x_iter[n], y_iter[n], ,]
    }

    # Optionally remove grids with all NA or empty
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
