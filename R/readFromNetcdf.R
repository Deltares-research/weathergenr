

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
readFromNetcdf <- function(
  in.path = NULL,
  in.file = NULL,
  dim.names = NULL,
  variables = NULL,
  origin.date = NULL)

  {

  ncin <- nc_open(paste0(in.path, in.file))

  # GET DIMENSIONS
  dim_names <- attributes(ncin$dim)$names
  DIMS <- vector(mode = "list", length = length(dim_names))
  names(DIMS) <- dim_names
  for (i in 1:length(dim_names)) {
    DIMS[[dim_names[i]]] <- ncvar_get(ncin,dim_names[i])
  }

  # GET VARIABLES
  VARS <- vector(mode = "list", length = length(variables))
  names(VARS) <- variables
  for (i in 1:length(variables)) {
    VARS[[variables[i]]] <- ncvar_get(ncin,variables[i])
  }

  # TRANSLATE INTO TIDY-FORMAT
  grid_mat <- VARS[[1]][,,1] %>% as_tibble() %>%
    setNames(1:length(DIMS[[dim.names$y]])) %>%
    mutate(xind = 1:n(), .before = 1) %>%
    gather(yind, data, -xind) %>%
    mutate(yind = as.numeric(yind)) %>%
    mutate(x = DIMS[[dim.names$x]][xind], y = DIMS[[dim.names$y]][yind], .after = yind) %>%
    na.omit() %>%
    mutate(data = list(NA))

  ##### Define variable data (as arrays)
  temp <- matrix(0, nrow = length(DIMS$time), ncol = length(variables)+1) %>%
    as_tibble() %>% setNames(c("date", variables)) %>%
    mutate(date = origin.date + DIMS[[dim.names$time]])

  for (n in 1:nrow(grid_mat)) {
    for (i in 1:length(variables)) {
      temp[,i+1] <- VARS[[i]][grid_mat$xind[n],grid_mat$yind[n], ]
    }
    grid_mat$data[[n]] <- temp
  }
  return(grid_mat)
}
