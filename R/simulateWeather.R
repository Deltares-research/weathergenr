
#' Main wrapper function for daily weather generation
#'
#' @param proj.name Placeholder.
#' @param wg.date.begin Placeholder.
#' @param wg.date.begin Placeholder.
#' @param wg.vars Placeholder.
#' @param wg.var.labs Placeholder.
#' @param wg.var.units Placeholder.
#' @param warm.var Placeholder.
#' @param ssig Placeholder.
#' @param ymax Placeholder.
#' @param rmax Placeholder.
#' @param nmax Placeholder.
#' @param nc.path Placeholder.
#' @param nc.file Placeholder.
#' @param nc.dimnames Placeholder.
#' @param validate Placeholder.
#' @param output.dir
#'
#' @return
#' @export
#' @import dplyr
#' @import tidyr
#' @import lubridate
#' @import tibble
#' @import patchwork
#' @importFrom sf st_as_sf st_sample st_cast st_coordinates
#' @importFrom readxl read_excel
simulateWeather <- function(
  proj.name = NULL,
  output.dir = NULL,
  climate_tidy = NULL,
  wg.date.begin = NULL,
  wg.vars = NULL,
  wg.var.labs = NULL,
  wg.var.units = NULL,
  warm.var = NULL,
  ssig = 0.90,
  ymax = 40,
  rmax = 20000,
  nmax = 5,
  nc.dimnames = NULL,
  validate = FALSE,
  ...)

 {

  # Number of grids
  grids  <- climate_tidy$id
  ngrids <- length(grids)

  #browser()
  cat("\u2059", "|", "Preprocessing historical data", "\r")

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

  cat("\u2713", "|", "Preprocessing historical data completed", "\r")


  #::::::::::: ANNUAL TIME-SERIES GENERATION USING WARM ::::::::::::::::::::::::

  message(cat("\u2059", "|", "Wavelet analysis of annual series", "\n"))

  #####  Wavelet analysis on observed annual series
  warm_annual_org <- climate_annual_wg_aavg %>% pull({{warm.var}})
  warm_annual <- warm_annual_org

  ####  Power spectra of observed annual series
  warm_power <- waveletAnalysis(variable = warm_annual, variable.unit = "mm",
    signif.level = ssig, plot = TRUE, out.path = out_path)

  ##### WAVELET DECOMPOSITION OF HISTORICAL SERIES
  wavelet_comps <- waveletDecompose(variable = warm_annual,
        signif.periods = warm_power$signif_periods,
        signif.level = ssig, plot = TRUE, out.path = out_path)

  message(cat("\u2713", "|", "Wavelet analysis: low frequency components defined", "\r"))

  message(cat("\u2059", "|", "Stochastic annual series generation started", "\n"))

  ##### Stochastic simuluation of historical annual series
  sim_annual <- waveletAR(wavelet.comps = wavelet_comps,
        num.years = ymax, num.realizations = rmax)

  message(cat("\u2713", "|", "Stochastic annual series generated completed", "\r"))


  message(cat("\u2059", "|", "Initial subsetting of stochastic series started", "\n"))

  # Wavelet analysis on simulated series
  sim_power <- sapply(1:rmax, function(x)
    waveletAnalysis(sim_annual[, x], signif.level = ssig)$GWS)

  ### Choose which parameters to consider for filtering
  sim_annual_sub <- waveletARSubset(
       series.obs = warm_annual,
       series.sim = sim_annual,
       power.obs = warm_power$GWS,
       power.sim = sim_power,
       power.period = warm_power$GWS_period,
       power.signif = warm_power$GWS_signif,
       nmax = nmax,
       out.path = out_path,
       ...)

  message(cat("\u2713", "|", "Stochastic traces selected:", ncol(sim_annual_sub$subsetted),
    "out of", ncol(sim_annual_sub$sampled), "meeting selection criteria", "\r"))

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

  sim_dates <- tibble(date = seq.Date(wg.date.begin, wg.date.begin-1 + years(ymax), by = "day")) %>%
    filter(!(month(date)==2 & day(date)==29)) %>% pull(date)

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

  message(cat("\u2713", "|", "Spatial-temporal dissaggregation completed"))

  #:::::::::::::::::::::: PERFORMANCE EVALUATION OF DAILY SIMULATED DATA :::::::

  if(isTRUE(validate)) {

    sample_ngrid <- 20

    sampleGrids <- sf::st_as_sf(climate_tidy[,c("x","y")], coords = c("x","y")) %>%
      sf::st_sample(size = sample_ngrid, type = "regular") %>%
      sf::st_cast("POINT") %>% st_coordinates() %>%
      as_tibble() %>%
      left_join(climate_tidy[,c("x","y","id")], by = c("X"="x","Y"="y")) %>%
      pull(id)

    sim_daily_wg_sample <- lapply(1:nmax, function(x) sim_daily_wg[[x]][sampleGrids])

    dailyPerformance(daily.sim = sim_daily_wg_sample,
                        daily.obs = climate_obs_daily[sampleGrids],
                        out.path = out_path,
                        variables = wg.vars[c(1,3,4)],
                        variable.labels = wg.var.labs[c(1,3,4)],
                        variable.units = wg.var.units[c(1,3,4)],
                        nmax = nmax)

    message(cat("\u2713", "|", "Validation plots generated and saved"))
  }



  #::::::::::::::::::::: OUTPUT HISTORICAL REALIZATIONS TO FILE ::::::::::::::::


  out_path_hist <- paste0(out_path, "hist_rlz/")

  # Dimensions in the outputted variable (order matters!)
  dim_ord <- names(nc_data$nc_dimensions)

  ncout_dims <- list()
  ncout_dims[[nc.dimnames$time]] <- ncdim_def(name = nc.dimnames$time,
      units = paste0("days since ", format(round(as.POSIXct(wg.date.begin), units = "day"),
        '%Y-%m-%d %M:%H:%S')), vals = 1:length(sim_dates), calendar = "no leap")

  ncout_dims[[nc.dimnames$y]] <- ncdim_def(nc.dimnames$y,
        units= "", vals = nc_data$nc_dimensions[[nc.dimnames$y]])

  ncout_dims[[nc.dimnames$x]] <- ncdim_def(name = nc.dimnames$x,
        units= "", vals = nc_data$nc_dimensions[[nc.dimnames$x]])

  ncout_dim_reorder <- unname(sapply(names(nc.dimnames), function(x) which(x == dim_ord)))

  # Other nc attributes
  ncout_varnames <- c(wg.vars, "pet")
  ncout_varunits <- c(wg.var.units, "mm/day")
  ncout_compression <- 4
  ncout_chunksize <- c(1,length(nc_data$nc_dimensions[[2]]),length(nc_data$nc_dimensions[[3]]))

  # ncout variables
  ncout_vars <- lapply(1:length(ncout_varnames), function(x) ncvar_def(
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
  names(ncout_vars) <- c(ncout_varnames, "spatial_ref")

  # TRANSLATE INTO TIDY-FORMAT
  coordGrid <- climate_tidy %>% mutate(date = list(NA))

  # template to store data from wg variables
  var_empty <- array(NA, c(
    ncout_dims[[nc.dimnames$time]]$len,
    ncout_dims[[nc.dimnames$y]]$len,
    ncout_dims[[nc.dimnames$x]]$len))

  #Loop through each variable and write data to netcdf
  ncout_vardata <- var_empty

  # Loop through each realization, calculate PET and save to file
  for (n in 1:nmax) {

    # create netCDF file and put arrays
    ncout_file <- nc_create(paste0(out_path_hist, proj.name,"_hist_rlz_",n, ".nc"),
      ncout_vars, force_v4 = TRUE)

    # New object to store the results
    daily_rlz <- sim_daily_wg[[n]]

    # Loop through each grid cell
    for (x in 1:ngrids) {

      # Calculate PET from temp, temp_min, temp_max
      daily_rlz[[x]]$pet <- with(daily_rlz[[x]],
          hargreavesPET(months = month(date), temp = temp,
            tdiff = temp_max - temp_min, lat = climate_tidy$y[x]))
    }

    #Loop through each variable and write data to netcdf
    for (i in 1:length(ncout_varnames)) {

      ncout_vardata[1:length(ncout_vardata)] <- var_empty

      for (c in 1:nrow(coordGrid)) {
        ncout_vardata[, coordGrid$yind[c], coordGrid$xind[c]] <- daily_rlz[[c]][[ncout_varnames[i]]]
      }

      # Put variables
      ncvar_put(ncout_file, varid = ncout_varnames[i], vals = ncout_vardata)

    }

    # Put spatial_def variable and attributes
    ncvar_put(ncout_file, varid = nc_data$nc_variables$spatial_ref$name,
              vals = nc_data$nc_variable_data$spatial_ref)

    sapply(1:length(nc_data$nc_attributes$spatial_ref), function(k)
      ncatt_put(ncout_file,
                varid = nc_data$nc_variables$spatial_ref$name,
                attname = names(nc_data$nc_attributes$spatial_ref)[k],
                attval = nc_data$nc_attributes$spatial_ref[[k]]))


   # Put global attributes
    if(length(nc_data$nc_attributes$global)>0) {

      sapply(1:length(nc_data$nc_attributes$global),
        function(k) ncatt_put(ncout_file, varid = 0,
                attname = names(nc_data$nc_attributes$global)[k],
                attval = nc_data$nc_attributes$global[[k]]))
    }

    nc_close(ncout_file)

  }
   message(cat("\u2713", "|", nmax, " stochastic realizations saved"))

}
