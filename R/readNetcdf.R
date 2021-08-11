

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
  origin.date = NULL)
  {


  get_dim_names <- function(f, v) {
    return(unlist(lapply(f$var[[v]]$dim, function(x) return(x$name))))
    }

  # Read-in file with netcdf4
  ncin <- nc_open(paste0(in.path, in.file))

  # GET DIMENSIONS
  dim_names <- attributes(ncin$dim)$names
  DIMS <- lapply(1:length(dim_names), function(x) ncvar_get(ncin,dim_names[x]))
  names(DIMS) <- dim_names

  # GET VARIABLES AND REORDER DIMENSIONS (x=lon,y=lat,time)
  VARS <- lapply(1:length(variables), function(x) ncvar_get(ncin,variables[x]))
  dim_ord <- sapply(variables, function(y) sapply(dim.names, function(x) which(x == get_dim_names(ncin, y))))
  VARS <- lapply(1:length(VARS), function(x) aperm(VARS[[x]], dim_ord[,x]))



  # translate variables to tidy format
  grid_mat <- VARS[[1]][,,1] %>%
    as_tibble(.name_repair = ~as.character(1:length(DIMS[[dim.names$y]]))) %>%
    mutate(xind = 1:n(), .before = 1) %>%
    gather(yind, data, -xind) %>%
    na.omit() %>%
    mutate(yind = as.numeric(yind)) %>%
    mutate(id = 1:n(), .before = xind) %>%
    mutate(x = DIMS[[dim.names$x]][xind],
           y = DIMS[[dim.names$y]][yind], .after = yind,
           data = list(NA))

  ##### Define variable data (as arrays)
  temp <- matrix(0, nrow = length(DIMS$time), ncol = length(variables)+1) %>%
    as_tibble(.name_repair = ~c("date", variables)) %>%
    mutate(date = origin.date + DIMS[[dim.names$time]])

  for (n in 1:nrow(grid_mat)) {
    temp[,-1] <- sapply(VARS, function(x) x[grid_mat$xind[n], grid_mat$yind[n], ])
    grid_mat$data[[n]] <- temp
  }
  return(grid_mat)
}
