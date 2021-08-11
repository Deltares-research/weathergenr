
#' Main wrapper function for daily weather generation
#'
#' @param proj.name Placeholder.
#' @param out_path Placeholder.
#' @param wg.date.begin Placeholder.
#' @param wg.date.end Placeholder.
#' @param sim.date.begin Placeholder.
#' @param wg.vars Placeholder.
#' @param wg.var.labs Placeholder.
#' @param wg.var.units Placeholder.
#' @param warm.var Placeholder.
#' @param ssig Placeholder.
#' @param ymax Placeholder.
#' @param rmax Placeholder.
#' @param nmax Placeholder.
#' @param PARCC Placeholder.
#' @param nc.path Placeholder.
#' @param nc.file Placeholder.
#' @param nc.dimnames Placeholder.
#' @param validate Placeholder.
#'
#' @return
#' @export
#' @import dplyr
#' @import tidyr
#' @import lubridate
#' @import readr
#' @import stringr
#' @import tibble
#' @import readxl
#' @import sf
simulateGriddedWeather <- function(
  proj.name = NULL,
  output.dir = NULL,
  wg.date.begin = NULL,
  sim.date.begin = NULL,
  wg.vars = NULL,
  wg.var.labs = NULL,
  wg.var.units = NULL,
  warm.var = NULL,
  ssig = 0.90,
  ymax = 40,
  rmax = 2000,
  nmax = 5,
  PARCC = NULL,
  nc.path = NULL,
  nc.file = NULL,
  nc.dimnames = NULL,
  perturb.file = NULL,
  validate = FALSE)

 {

  # ::::::::: Load & prepare climate data ::::::::::::::::::::::::::::::::::::::

  # GABON/ ERA5-GRIDDED CLIMATE DATA
  climate_tidy <- readNetcdf(
    in.path = nc.path,
    in.file = nc.file,
    dim.names = list(x = "lon", y = "lat", time = "time"),
    variables = wg.vars,
    origin.date = wg.date.begin)

  #Remove extra grids
  grid_select <- which(sapply(1:nrow(climate_tidy), function(x)
    !is.na(mean(climate_tidy$data[[x]]$temp_min))))

  climate_tidy <- climate_tidy[grid_select,] %>%
    mutate(id = 1:n())

  # Number of grids
  grids  <- climate_tidy$id;
  ngrids <- length(grids)

  message(cat("\u2713", "|", "Loaded historical climate: ", ngrids, "grids"))

  ##### Set historical gridded data to be inputted to wg (remove leap days!)
  noleap_index <- with(climate_tidy$data[[1]], which(!(month(date) == 2 & day(date) == 29)))
  climate_obs_daily <- lapply(climate_tidy$data, function(x) x[noleap_index, ])
  wg_dates <- climate_obs_daily[[1]]$date

  # Month order for water year
  wy_months <- setdiff(c(month(wg.date.begin):12, 1:(month(wg.date.begin)-1)),0)

  # Date template table with water-year adjustment
  wg_dates_adjusted <- tibble(year = year(wg_dates) - min(year(wg_dates)) + 1,
      mon = month(wg_dates), day = day(wg_dates)) %>%
      mutate(mon = as.numeric(as.character(factor(mon, levels = wy_months, labels= 1:12))))

  # Date vectors for the final simulated series
  sim_dates_leap <- seq.Date(sim.date.begin, sim.date.begin + years(ymax)-1, by = "day")
  sim_dates <- sim_dates_leap[!(month(sim_dates_leap)==2 & day(sim_dates_leap)==29)]

  # Prepare tables for daily, monthly, and annual values
  climate_daily_wg <- lapply(1:ngrids, function(i)
        bind_cols(wg_dates_adjusted, climate_obs_daily[[i]][,-1]) %>%
        arrange(year, mon, day))

  climate_monthly_wg <- lapply(1:ngrids, function(i)
        climate_daily_wg[[i]] %>% group_by(year, mon) %>%
        summarize(across({{wg.vars}}, mean)) %>%
        ungroup() %>% suppressMessages())

  climate_annual_wg <- lapply(1:ngrids, function(i)
        climate_monthly_wg[[i]] %>%
        group_by(year) %>%
        summarize(across({{wg.vars}}, mean)) %>%
        ungroup() %>% suppressMessages())

  # Obtain scalar (area-averaged) weather time-series for each variable
  climate_daily_wg_aavg   <- Reduce(`+`, climate_daily_wg) / ngrids
  climate_monthly_wg_aavg <- Reduce(`+`, climate_monthly_wg) / ngrids
  climate_annual_wg_aavg  <- Reduce(`+`, climate_annual_wg) / ngrids

  # Matrix of monthly historical means (rows = years, columns=month 1 to 12)
  wg_variable_matrix <- climate_monthly_wg_aavg %>%
    dplyr::select(year, mon, {{warm.var}}) %>%
    spread(mon, {{warm.var}}) %>%
    dplyr::select(-year)

  # Matrix of monthly historical means normalized over the period
  wg_monthly_factors <- t(apply(wg_variable_matrix, 1, function(x) x/mean(x)))

  # Indices of all days (rows = year & mon combination, columns= days 1 to 31)
  wg_day_of_year <- climate_daily_wg_aavg %>%
    dplyr::select(year, mon, day) %>%
    mutate(ind = 1:n()) %>%
    spread(key = day, value = ind)

  #::::::::::: ANNUAL TIME-SERIES GENERATION USING WARM ::::::::::::::::::::::::

  #####  Wavelet analysis on observed annual series
  warm_annual_org <- climate_annual_wg_aavg %>% pull({{warm.var}})
  warm_annual <- warm_annual_org

  ####  Power spectra of observed annual series
  warm_power <- waveletAnalysis(variable = warm_annual,
         variable.unit = "mm", signif.level = ssig)

  # Save results to file
  ggsave(plot = warm_power$plot, height = 9, width = 12,
         filename = paste0(out_path, "warm_observed_wavelet.png"))

  message(cat("\u2713", "|", "Wavelet Analysis"))

  ##### WAVELET DECOMPOSITION OF HISTORICAL SERIES
  wavelet_comps <- waveletDecompose(variable = warm_annual,
        signif.periods = warm_power$signif_periods,
        signif.level = ssig)

  message(cat("\u2713", "|", "Wavelet Decomposition"))

  ##### Stochastic simuluation of historical annual series
  sim_annual <- waveletAR(wavelet.comps = wavelet_comps$out,
                           num.years = ymax,
                           num.realizations = rmax)

  message(cat("\u2713", "|", "Wavelet AR Model simulation"))

  # Wavelet analysis on simulated series
  sim_power <- sapply(1:rmax, function(x)
    waveletAnalysis(sim_annual[, x], signif.level = ssig)$GWS)

  #p <- waveletPlot(power.period = warm_power$GWS_period,
  #    power.signif = warm_power$GWS_signif, power.obs = warm_power$GWS,
  #    power.sim = sim_power)

  ### Choose which parameters to consider for filtering
  sim_annual_sub <- waveletARSubset(
       series.obs = warm_annual,
       series.sim = sim_annual,
       power.obs = warm_power$GWS,
       power.sim = sim_power,
       power.period = warm_power$GWS_period,
       power.signif = warm_power$GWS_signif,
       mean.bounds = c(0.97, 1.03),
       sdev.bounds = c(0.98, 1.05),
       max.bounds  = c(0.95, 1.10),
       min.bounds  = c(0.90, 1.05),
       power.bounds = c(0.70, 2),
       nonsig.threshold = 0.8,
       nmax = nmax,
       save.plots = TRUE,
       save.series = TRUE,
       verbose = TRUE,
       out.path = out_path)

  #::::TEMPORAL/SPATIAL DISSAGGREGATION ::::::::::::::::::::::::::::::::::::::::::

  # Subsetted realizations of annual simulated time-series
  sim_annual_final <- sim_annual_sub$sampled

  # For each simulated annual value, find the year with closest value
  closest_yr_index <- apply(sim_annual_final, c(1,2), function(x)
    which.min(abs(x - warm_annual)))

  #Monthly factors
  sim_monthly_ini <- lapply(1:nmax, function(x)
    wg_monthly_factors[closest_yr_index[,x],] * sim_annual_final[,x])

  #:::::::::::::::::::::: KNN-SAMPLING :::::::::::::::::::::::::::::::::::::::::

  sim_daily_wg <- vector(mode = "list", length = nmax)
  knn_optimk <- ceiling(sqrt(length(warm_annual)))

  for (n in 1:nmax) {

    #Obtain KNN table with sampled years for each month
    knn_sample <- sapply(1:12, function(m) sapply(1:ymax, function(y)
      kNearestNeighbors(sim_monthly_ini[[n]][y,m], knn_optimk,
                        wg_variable_matrix[,m]))) %>%
        as.data.frame() %>% setNames(1:12) %>%
        mutate(syear = 1:n()) %>%
      gather(key = mon, value = year, -syear) %>%
      mutate(mon = as.numeric(mon)) %>%
      arrange(syear, mon)

    # Find the indices of sampled days
    rlz_index <- knn_sample %>%
      left_join(wg_day_of_year, by = c("year","mon")) %>%
      gather(key = day, value = value, -syear, -mon, -year) %>%
      arrange(syear, mon) %>% dplyr::select(-year) %>% na.omit() %>% pull(value)

    rlz_index2 <- knn_sample %>%
      left_join(wg_day_of_year, by = c("year","mon")) %>%
      gather(key = day, value = value, -syear, -mon, -year) %>%
      arrange(syear, mon) %>% dplyr::select(-year) %>% pull(value)


    # Construct final dataset
    sim_daily_wg[[n]] <- lapply(1:ngrids, function(x) {
        climate_daily_wg[[x]][rlz_index, ] %>%
            dplyr::select({{wg.vars}}) %>%
            mutate(date = sim_dates, .before = 1)
    })

  }

  message(cat("\u2713", "|", "Spatial-temporal dissaggregation of annual simulated series"))

  #:::::::::::::::::::::: PERFORMANCE EVALUATION OF DAILY SIMULATED DATA :::::::

  if(isTRUE(validate)) {

    sampleGrids <- sf::st_as_sf(climate_tidy[,c("x","y")], coords = c("x","y")) %>%
      st_sample(size=20, type = "regular") %>%
      st_cast("POINT") %>% st_coordinates() %>%
      as_tibble() %>%
      left_join(climate_tidy[,c("x","y","id")], by = c("X"="x","Y"="y")) %>%
      pull(id)

    sim_daily_wg_sample <- lapply(1:nmax, function(x) sim_daily_wg[[x]][sampleGrids])

    dailyPerformance(daily.sim = sim_daily_wg_sample,
                        daily.obs = climate_obs_daily[sampleGrids],
                        out.path = out_path,
                        variables = wg.vars[c(1,3,4)],
                        variable.labels = wg.var.labs[c(1,3,4)],
                        variable.units = wg.var.units[c(1,3,4)])
  }


  #:::::::::::::::::::::: CLIMATE CHANGE PERTURBATIONS :::::::::::::::::::::::::::

  perturb1 <- read_excel(perturb.file, sheet = "par1", range = "B3:D36")
  perturb2 <- read_excel(perturb.file, sheet = "par2", range = "B3:D36")

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
  scn_mat <- expand_grid(nvar=1:nmax, par1=1:PARCC[[1]]$increments, par2=1:PARCC[[2]]$increments)
  smax <- nrow(scn_mat)

  message(cat("\u2713", "|", "Climate change scenario matrix created"))

  #::::::::::::::::::::::: TEMPLATE FOR WRITING TO NETCDF ::::::::::::::::::::::

  nc_varnames <- c(wg.vars, "pet")
  ncvar_units <- c(wg.var.units, "mm/day")
  compression <- 4

  # Open template netcdf file
  nc_in <- nc_open(paste0(nc.path, nc.file))

  # Global attributes
  nc_atts <- ncatt_get(nc_in, 0)

   # Put spatial variable attributes
  spatial_ref_name <- "spatial_ref"
  spatial_ref_val <- ncvar_get(nc_in, spatial_ref_name)
  spatial_ref_attribs <- ncatt_get(nc_in,spatial_ref_name)
  spatial_ref_def <- ncvar_def(name = spatial_ref_name, units = nc_in$var[[spatial_ref_name]]$units,
    dim = nc_in$var[[spatial_ref_name]]$dims, prec = "integer")

  #::::::::::::::::::::::: NC DIMENSIONS :::::::::::::::::::::::::::::::::::::::

  nc_dims <- list()
  nc_dims[[nc.dimnames$x]] <- ncdim_def(nc.dimnames$x,
        units= "", vals = ncvar_get(nc_in,nc.dimnames$x))

  nc_dims[[nc.dimnames$y]] <- ncdim_def(nc.dimnames$y,
        units= "", vals = ncvar_get(nc_in,nc.dimnames$y))

  # Set time dimension
  nc_dims[[nc.dimnames$time]] <- ncdim_def(nc.dimnames$time,
      units = paste0("days since ", format(round(as.POSIXct(sim.date.begin), units = "day"),
        '%Y-%m-%d %M:%H:%S')),
      vals = 1:length(sim_dates),
      calendar = "noleap")
  names(nc_dims) <- unlist(nc.dimnames, use.names = FALSE)


  #::::::::::::::::::::: NC VARIABLES ::::::::::::::::::::::::::::::::::::::::::

  nc_vars <- lapply(1:length(nc_varnames), function(x) ncvar_def(
        name = nc_varnames[x],
        units = ncvar_units[x],
        dim = nc_dims,
        missval = NA,
        longname = nc_varnames[x],
        prec =  "float",
        shuffle = FALSE,
        compression = compression,
        chunksizes=c(nc_dims[[1]]$len,nc_dims[[2]]$len,1)))

  nc_vars_all <- nc_vars
  nc_vars_all[[length(nc_vars_all)+1]] <- spatial_ref_def

  # TRANSLATE INTO TIDY-FORMAT
  coordGrid <- climate_tidy %>% mutate(date = list(NA))

  # template to store data from wg variables
  var_empty <- array(NA, c(nc_dims[[nc.dimnames$x]]$len,
        nc_dims[[nc.dimnames$y]]$len, nc_dims[[nc.dimnames$time]]$len))

  #::::::::::::::::::: LOOP THROUGH CLIMATE CHANGES ::::::::::::::::::::::::::::

  # Loop through each scenario
  for (s in 1:smax) {

    # New object to store the results
    daily_rlz <- sim_daily_wg[[scn_mat$nvar[s]]]

    # Current perturbation scenario for each variable
    perturb_par1 <- list(mean = PARCC[[1]]$mean$steps[scn_mat$par1[s],],
           var = PARCC[[1]]$var$steps[scn_mat$par1[s],])

    perturb_par2 <- list(mean = PARCC[[2]]$mean$steps[scn_mat$par2[s],],
           var = PARCC[[2]]$var$steps[scn_mat$par2[s],])

    # Loop through each grid cell
    for (x in 1:ngrids) {

      daily_rlz[[x]]$precip <- quantileMapping(value = daily_rlz[[x]]$precip,
            date = sim_dates, par = perturb_par1, operator = "multiply")

      daily_rlz[[x]]$temp <- quantileMapping(value = daily_rlz[[x]]$temp,
            date = sim_dates, par = perturb_par2, operator = "add")

      daily_rlz[[x]]$temp_min <- quantileMapping(value = daily_rlz[[x]]$temp_min,
            date = sim_dates, par = perturb_par2, operator = "add")

      daily_rlz[[x]]$temp_max <- quantileMapping(value = daily_rlz[[x]]$temp_max,
            date = sim_dates, par = perturb_par2, operator = "add")

      daily_rlz[[x]]$pet <- with(daily_rlz[[x]],
          hargreavesPET(months = month(date), temp = temp,
            tdiff = temp_max - temp_min, lat = climate_tidy$y[x]))
    }

    # create netCDF file and put arrays
    nc_out <- nc_create(paste0(output.dir, proj.name,"_rlz_",s, ".nc"), nc_vars_all, force_v4 = TRUE)

    # Put spatial_def variable
    ncvar_put(nc_out, varid = spatial_ref_def, vals = spatial_ref_val)

      # Put variable attributes
    sapply(1:length(spatial_ref_attribs), function(k) ncatt_put(nc_out,
      varid = spatial_ref_name, attname = names(spatial_ref_attribs)[k],
      attval = spatial_ref_attribs[[k]]))

   # Put global attributes
    if(length(nc_atts)>0) {
      sapply(1:length(nc_atts), function(k)
        ncatt_put(nc_out, varid = 0, attname = names(nc_atts)[k],
                  attval = nc_atts[[k]]))
    }

    #Loop through each variable and write data to netcdf
    nc_var_data <- var_empty

    for (i in 1:length(nc_vars)) {

      nc_var_data[1:length(nc_var_data)] <- var_empty

      for (c in 1:nrow(coordGrid)) {
        nc_var_data[coordGrid$xind[c],coordGrid$yind[c], ] <- daily_rlz[[c]][[nc_varnames[i]]]
      }

      # Put variables
      ncvar_put(nc_out, varid = nc_vars[[i]], vals = nc_var_data)

    }
    nc_close(nc_out)
  }

  message(cat("\u2713", "|", "Results outputted to", s, "netcdf files at: ", out_path))
}


