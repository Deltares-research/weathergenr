

  ##############################################################################
  ##############################################################################


#' Climate Change perturbations function
#'
#' @param proj.name To be completed...
#' @param input.dir  To be completed...
#' @param out.path To be completed...
#' @param sim.date.begin To be completed...
#' @param wg.vars To be completed...
#' @param wg.var.labs To be completed...
#' @param wg.var.units To be completed...
#' @param nc.dimnames To be completed...
#' @param change.settings To be completed...
#'
#' @return
#' @export
#' @import dplyr
#' @import ncdf4
#' @import lubridate
imposeClimateChanges <- function(
  proj.name = NULL,
  in.path = NULL,
  in.file = NULL,
  file.suffix = "",
  out.path = NULL,
  sim.date.begin = NULL,
  wg.vars = NULL,
  wg.var.units = NULL,
  nc.dimnames = NULL,
  change.settings = NULL,
  save.scenario.matrix = TRUE)

{

  #:::::::::::::::::::::: CLIMATE CHANGE SETTINGS ::::::::::::::::::::::::::::::

  perturb1 <- readxl::read_excel(change.settings, sheet = "par1", range = "B3:D36")
  perturb2 <- readxl::read_excel(change.settings, sheet = "par2", range = "B3:D36")

  PARCC <- list(par1 = list(), par2 = list())

  PARCC$par1 <- list(
    name = perturb1[[1,2]],
    change_op <- perturb1[[2,2]],
    increments = as.numeric(perturb1[[3,2]]),
    mean = list(min = as.numeric(perturb1[7:18,1:3] %>% pull(2)),
                max = as.numeric(perturb1[7:18,1:3] %>% pull(3)),
                obs = as.numeric(perturb1[7:18,1:3] %>% pull(1))),
    var  = list(min = as.numeric(perturb1[22:33,1:3] %>% pull(2)),
                max = as.numeric(perturb1[22:33,1:3] %>% pull(3)),
                obs = as.numeric(perturb1[22:33,1:3] %>% pull(1)))
  )
  PARCC$par2 <- list(
    name = perturb2[[1,2]],
    change_op <- perturb2[[2,2]],
    increments = as.numeric(perturb2[[3,2]]),
    mean = list(min = as.numeric(perturb2[7:18,1:3] %>% pull(2)),
                max = as.numeric(perturb2[7:18,1:3] %>% pull(3)),
                obs = as.numeric(perturb2[7:18,1:3] %>% pull(1))),
    var  = list(min = as.numeric(perturb2[22:33,1:3] %>% pull(2)),
                max = as.numeric(perturb2[22:33,1:3] %>% pull(3)),
                obs = as.numeric(perturb2[22:33,1:3] %>% pull(1)))
  )

  PARCC$par1$sind <- 1:PARCC$par1$increments
  PARCC$par2$sind <- 1:PARCC$par2$increments


  for (i in 1:length(PARCC)) {

    PARCC[[i]]$mean$steps <- sapply(1:12, function(m)
      seq(PARCC[[i]]$mean$min[m], PARCC[[i]]$mean$max[m], length.out = PARCC[[i]]$increments)) %>%
      rbind(rep(NA,12))

    PARCC[[i]]$var$steps <-   sapply(1:12, function(m)
      seq(PARCC[[i]]$var$min[m], PARCC[[i]]$var$max[m], length.out = PARCC[[i]]$increments)) %>%
      rbind(rep(NA,12))
  }

  # Scenario matrix
  scn_mat <- expand_grid(par1 = 1:PARCC[[1]]$increments, par2 = 1:PARCC[[2]]$increments) %>%
    mutate(id = 1:n(), .before = 1) %>%
    mutate(precip_change_mean = PARCC[[1]]$mean$steps[par1,1], temp_change_mean = PARCC[[2]]$mean$steps[par2,1]) %>%
    mutate(precip_change_variance = PARCC[[1]]$var$steps[par1,1], temp_change_variance = PARCC[[2]]$var$steps[par2,1])

  if(isTRUE(save.scenario.matrix)) {
    write.csv(x = scn_mat, file = paste0(out.path, "scenario_matrix.csv"))
  }

  smax <- nrow(scn_mat)

  message(cat("\u2713", "|", "scenario matrix created: ", smax, "scenarios in total"))

  #::::::::::::::::::::::: TEMPLATE FOR WRITING TO NETCDF ::::::::::::::::::::::

  nc_data <- readNetcdf(
    in.path = in.path,
    in.file = in.file,
    dim.names = list(x = "lon", y = "lat", time = "time"),
    variables = wg.vars,
    origin.date = sim.date.begin,
    leap.year = FALSE)

  # TRANSLATE INTO TIDY-FORMAT
  coordGrid <- nc_data$tidy_data %>% dplyr::select(-data)

  # Number of grids
  grids  <- coordGrid$id
  ngrids <- length(grids)


   sim_dates <- nc_data$tidy_data$data[[1]]$date

  # Dimensions in the outputted variable (order matters!)
  dim_ord <- names(nc_data$nc_dimensions)

  ncout_dims <- list()
  ncout_dims[[nc.dimnames$time]] <- ncdf4::ncdim_def(name = nc.dimnames$time,
      units = paste0("days since ", format(round(as.POSIXct(sim.date.begin), units = "day"),
        '%Y-%m-%d %M:%H:%S')), vals = 1:length(sim_dates), calendar = "no leap")

  ncout_dims[[nc.dimnames$y]] <- ncdf4::ncdim_def(nc.dimnames$y,
        units= "", vals = nc_data$nc_dimensions[[nc.dimnames$y]])

  ncout_dims[[nc.dimnames$x]] <- ncdf4::ncdim_def(name = nc.dimnames$x,
        units= "", vals = nc_data$nc_dimensions[[nc.dimnames$x]])

  ncout_dim_reorder <- unname(sapply(names(nc.dimnames), function(x) which(x == dim_ord)))

  # Other nc attributes
  ncout_varnames <- c(wg.vars, "pet")
  ncout_varunits <- c(wg.var.units, "mm/day")
  ncout_compression <- 4
  ncout_chunksize <- c(1,length(nc_data$nc_dimensions[[2]]),length(nc_data$nc_dimensions[[3]]))

  # ncout variables
  ncout_vars <- lapply(1:length(ncout_varnames), function(x) ncdf4::ncvar_def(
        name = ncout_varnames[x],
        units = ncout_varunits[x],
        dim = ncout_dims,
        missval = NA,
        longname = ncout_varnames[x],
        prec =  "float",
        shuffle = FALSE,
        compression = ncout_compression,
        chunksizes= ncout_chunksize))

  ncout_vars[[length(ncout_vars)+1]] <- nc_data$nc_variables$spatial_ref


  # template to store data from wg variables
  var_empty <- array(NA, c(
    ncout_dims[[nc.dimnames$time]]$len,
    ncout_dims[[nc.dimnames$y]]$len,
    ncout_dims[[nc.dimnames$x]]$len))

  #Loop through each variable and write data to netcdf
  ncout_vardata <- var_empty



  #::::::::::::::::::: LOOP THROUGH CLIMATE CHANGES ::::::::::::::::::::::::::::

  #Create output directory if doesn't exist
  if (!dir.exists(out.path)) {dir.create(out.path)}

  # Loop through each scenario
  for (s in 1:smax) {

    ncout_file <- ncdf4::nc_create(paste0(out.path, proj.name,"_climate_change_",file.suffix,"_", s, ".nc"),
      ncout_vars, force_v4 = TRUE)

    # New object to store the results
    daily_rlz <- nc_data$tidy_data$data

    # Current perturbation scenario for each variable
    perturb_par1 <- list(mean = PARCC[[1]]$mean$steps[scn_mat$par1[s],],
           var = PARCC[[1]]$var$steps[scn_mat$par1[s],])

    perturb_par2 <- list(mean = PARCC[[2]]$mean$steps[scn_mat$par2[s],],
           var = PARCC[[2]]$var$steps[scn_mat$par2[s],])


    # Loop through each grid cell
    for (x in 1:ngrids) {

      # Perturb daily precipitation
      daily_rlz[[x]]$precip <- quantileMapping(value = daily_rlz[[x]]$precip,
            date = sim_dates, par = perturb_par1, operator = "multiply")


      #################################################################################
      # Perturb temp, temp_min, and temp_max
      daily_rlz[[x]]$temp <- quantileMapping(value = daily_rlz[[x]]$temp,
            date = sim_dates, par = perturb_par2, operator = "add")
      daily_rlz[[x]]$temp_min <- quantileMapping(value = daily_rlz[[x]]$temp_min,
            date = sim_dates, par = perturb_par2, operator = "add")
      daily_rlz[[x]]$temp_max <- quantileMapping(value = daily_rlz[[x]]$temp_max,
            date = sim_dates, par = perturb_par2, operator = "add")

      #################################################################################

      # Calculate PET from temp, temp_min, temp_max
      daily_rlz[[x]]$pet <- with(daily_rlz[[x]],
          hargreavesPET(months = month(date), temp = temp,
            tdiff = temp_max - temp_min, lat = coordGrid$y[x]))
    }

    #Loop through each variable and write data to netcdf
    for (i in 1:length(ncout_varnames)) {

      ncout_vardata <- var_empty

      for (x in 1:ngrids) {
        ncout_vardata[, coordGrid$yind[x], coordGrid$xind[x]] <- daily_rlz[[x]][[ncout_varnames[i]]]
      }

      # Put variables
      ncdf4::ncvar_put(ncout_file, varid = ncout_varnames[i], vals = ncout_vardata)
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

  message(cat("\u2713", "|", "Results outputted to", s, "netcdf files at: ", out_path))

}

