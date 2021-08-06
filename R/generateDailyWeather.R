
#' Main wrapper function for daily weather generation
#'
#' @param project_name Placeholder.
#' @param out_path Placeholder.
#' @param wg_date_begin Placeholder.
#' @param wg_date_end Placeholder.
#' @param sim_date_begin Placeholder.
#' @param wg_variables Placeholder.
#' @param wg_variable_labs Placeholder.
#' @param wg_variable_units Placeholder.
#' @param warm_variable Placeholder.
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
#'
generateDailyWeather <- function(
  project_name = NULL,
  wg_date_begin = NULL,
  wg_date_end = NULL,
  sim_date_begin = NULL,
  wg_variables = NULL,
  wg_variable_labs = NULL,
  wg_variable_units = NULL,
  warm_variable = NULL,
  ssig = 0.90,
  ymax = 40,
  rmax = 10000,
  nmax = 5,
  PARCC = NULL,
  nc.path = NULL,
  nc.file = NULL,
  nc.dimnames = NULL,
  validate = FALSE
  )

 {

  out_path = paste0("./areas/", project_name,"/output/")

  # GABON/ ERA5-GRIDDED CLIMATE DATA
  climate_input_tidy <- readFromNetcdf(
    in.path = nc.path,
    in.file = nc.file,
    variables = wg_variables,
    dim.names = nc.dimnames,
    origin.date = origin.date)


  # Number of grid cells (or point-locations)
  grids  <- climate_input_tidy$id
  ngrids <- length(grids)

  message(cat("\u2713", "|", "Loaded historical climate data: total grids:", ngrids))

  # READ-IN & PROCESS CLIMATE DATA -----------------------------------------------

  # Create main project dir (if not already existing)
  ifelse(!dir.exists(file.path(paste0("./areas/", project_name))),
         dir.create(file.path(paste0("./areas/", project_name))), FALSE)

  # Create input directory (if not already existing)
  ifelse(!dir.exists(file.path(paste0("./areas/", project_name,"./input/"))),
         dir.create(file.path(paste0("./areas/", project_name,"./input/"))), FALSE)

  # Create output directory directory (if not already existing)
  ifelse(!dir.exists(file.path(paste0("./areas/", project_name,"./output/"))),
         dir.create(file.path(paste0("./areas/", project_name,"./output/"))), FALSE)

  climate_input_tidyX <- climate_input_tidy
  climate_input_tidyX$data <- lapply(1:nrow(climate_input_tidy),
    function(x) filter(climate_input_tidy$data[[x]], date %in% wg_obs_dates))


  ##### Set historical gridded data to be inputted to wg
  climate_obs_daily <- climate_input_tidyX$data

  #wg_date_begin <- as.Date(climate_input_tidy$data[[1]]$date[1])
  #wg_date_end <- as.Date(climate_input_tidy$data[[1]]$date[length(climate_input_tidy$data[[1]]$date)])

  # Date vectors for the historical input series
  wg_obs_dates_leap <- seq.Date(wg_date_begin, wg_date_end, by = "day") #temp var
  wg_obs_dates <- wg_obs_dates_leap[!(month(wg_obs_dates_leap)==2 & day(wg_obs_dates_leap)==29)]
  wg_obs_year_num <- length(wg_obs_dates %>% floor_date("month") %>% unique())/12

  # Date vectors for the final simulated series
  wg_sim_dates_leap <- seq.Date(sim_date_begin, sim_date_begin + years(ymax)-1, by = "day")
  wg_sim_dates <- wg_sim_dates_leap[!(month(wg_sim_dates_leap)==2 & day(wg_sim_dates_leap)==29)]

  # Month order for water year
  wg_month_order <- setdiff(c(month(wg_date_begin):12, 1:(month(wg_date_begin)-1)),0)

  ############## FIGURE OUT HOW TO START AT A RANDOM DAY?????????

  #### Eliminate leap days and calculate daily, monthly, and annual values for wg
  climate_obs_daily_wg <- vector("list", ngrids)
  climate_obs_monthly_wg <- vector("list", ngrids)
  climate_obs_annual_wg <- vector("list", ngrids)
  for (i in 1:ngrids) {

    climate_obs_daily_wg[[i]] <- climate_obs_daily[[i]] %>%
      filter(date %in% wg_obs_dates) %>%
      mutate(year = rep(1:wg_obs_year_num, each = 365),
             mon = as.numeric(as.character(factor(month(date),
                  levels = wg_month_order, labels= 1:12))),
             day = day(date)) %>%
      dplyr::select(year, mon, day, {{wg_variables}}) %>%
      arrange(year, mon, day)

    climate_obs_monthly_wg[[i]] <- climate_obs_daily_wg[[i]] %>%
      group_by(year, mon) %>%
      summarize(across({{wg_variables}}, mean)) %>%
      ungroup() %>% suppressMessages()

    climate_obs_annual_wg[[i]] <- climate_obs_daily_wg[[i]] %>%
      group_by(year) %>%
      summarize(across({{wg_variables}}, mean)) %>%
      ungroup() %>% suppressMessages()

  }

  # Obtain scalar (area-averaged) weather time-series for each variable
  climate_obs_daily_wg_aavg   <- Reduce(`+`, climate_obs_daily_wg) / ngrids
  climate_obs_monthly_wg_aavg <- Reduce(`+`, climate_obs_monthly_wg) / ngrids
  climate_obs_annual_wg_aavg  <- Reduce(`+`, climate_obs_annual_wg) / ngrids

  # Matrix of historical data (rows = years 1 to ymax, columns = months 1 to 12)
  wg_obs_mat <- climate_obs_monthly_wg_aavg %>%
    dplyr::select(year, mon, {{warm_variable}}) %>%
    spread(mon, {{warm_variable}}) %>%
    dplyr::select(-year)

  # Monthly factors for hist. data (rows = 1 to ymax, columns = 1 to 12)
  wg_obs_monthly_factors <- t(apply(wg_obs_mat, 1, function(x) x/mean(x)))

  # Indices of all days (rows = year & mon combination, columns= days 1 to 31)
  wg_dayofyear <- climate_obs_daily_wg_aavg %>%
    dplyr::select(year, mon, day) %>%
    mutate(ind = 1:n()) %>%
    spread(key = day, value = ind)

  #::::::::::: ANNUAL TIME-SERIES GENERATION USING WARM ::::::::::::::::::::::::::

  warm_obs_annual_org <- climate_obs_annual_wg_aavg %>% pull({{warm_variable}})

  #####  Wavelet analysis on observed annual series
  warm_obs_annual <- warm_obs_annual_org

  ####  Power spectra of observed annual series
  warm_obs_power <- waveletAnalysis(variable = warm_obs_annual,
                                     variable.unit = "mm",
                                     signif.level = ssig,
                                     save.plot = TRUE)

  message(cat("\u2713", "|", "Wavelet Analysis"))

  ##### WAVELET DECOMPOSITION OF HISTORICAL SERIES
  wavelet_comps <- waveletDecomposition(variable = warm_obs_annual,
                            signif.period.list = warm_obs_power$signif_periods_list,
                            signif.level = ssig)

  message(cat("\u2713", "|", "Wavelet Decomposition"))

  ##### Stochastic simuluation of historical annual series
  warm_sim_annual <- waveletAR(wavelet.comps = wavelet_comps$out,
                           num.years = ymax,
                           num.realizations = rmax)

  message(cat("\u2713", "|", "Wavelet AR Model simulation"))

  # Wavelet analysis on simulated series
  warm_sim_power <- sapply(1:rmax, function(x)
    waveletAnalysis(warm_sim_annual[, x], signif.level = ssig)$GWS)

  p <- waveletGlobalSpectrumPlot(power.period = warm_obs_power$GWS_period,
                 power.signif = warm_obs_power$GWS_signif,
                 power.obs = warm_obs_power$GWS,
                 power.sim = warm_sim_power)

    ### Choose which parameters to consider for filtering
  warm_sim_annual_sub <- waveletARSubset(
       series.obs = warm_obs_annual,
       series.sim = warm_sim_annual,
       power.obs = warm_obs_power$GWS,
       power.sim = warm_sim_power,
       power.period = warm_obs_power$GWS_period,
       power.signif = warm_obs_power$GWS_signif,
       mean.bounds = c(0.90, 1.15),
       sdev.bounds = c(0.85, 1.15),
       max.bounds  = c(0.85, 1.15),
       min.bounds  = c(0.80, 1.15),
       power.bounds = c(0.70, 5),
       nonsig.threshold = 2,
       nmax = nmax,
       save.plots = TRUE,
       save.series = TRUE,
       verbose = FALSE,
       out.path = out_path)

  #::::TEMPORAL/SPATIAL DISSAGGREGATION ::::::::::::::::::::::::::::::::::::::::::

  # Subsetted realizations of annual simulated time-series
  warm_sim_annual_final <- warm_sim_annual_sub$sampled

  # For each simulated annual value, find the year with closest value
  closest_yr_index <- apply(warm_sim_annual_final, c(1,2), function(x)
    which.min(abs(x - warm_obs_annual)))

  #Monthly factors
  warm_sim_monthly_ini <- lapply(1:nmax, function(x)
    wg_obs_monthly_factors[closest_yr_index[,x],] * warm_sim_annual_final[,x])

  #:::::::::::::::::::::: KNN-SAMPLING :::::::::::::::::::::::::::::::::::::::::::

  climate_sim_daily_wg <- vector(mode = "list", length = nmax)
  knn_optimk <- ceiling(sqrt(length(warm_obs_annual)))
  for (n in 1:nmax) {

    #Obtain KNN table with sampled years for each month
    knn_sample <- sapply(1:12, function(m) sapply(1:ymax, function(y)
      kNearestNeighbors(warm_sim_monthly_ini[[n]][y,m], knn_optimk,
                        wg_obs_mat[,m]))) %>%
        as.data.frame() %>% setNames(1:12) %>%
        mutate(syear = 1:n()) %>%
      gather(key = mon, value = year, -syear) %>%
      mutate(mon = as.numeric(mon)) %>%
      arrange(syear, mon)

    # Find the indices of sampled days
    rlz_index <- knn_sample %>%
      left_join(wg_dayofyear, by = c("year","mon")) %>%
      gather(key = day, value = value, -syear, -mon, -year) %>%
      arrange(syear, mon) %>% dplyr::select(-year) %>% na.omit() %>% pull(value)

    # Construct final dataset
    climate_sim_daily_wg[[n]] <- lapply(1:ngrids, function(x) {
        climate_obs_daily_wg[[x]][rlz_index, ] %>%
            dplyr::select({{wg_variables}}) %>%
            mutate(date = wg_sim_dates, .before = 1)
    })

  }

  message(cat("\u2713", "|", "Spatial-temporal dissaggregation of annual simulated series"))

  #:::::::::::::::::::::: PERFORMANCE EVALUATION OF DAILY SIMULATED DATA :::::::::

  if(isTRUE(validate)) {

    daily_obs_check <- (lapply(1:ngrids, function(x)
    climate_obs_daily[[x]] %>% filter(date %in% wg_obs_dates)))

    daily_sim_check <- lapply(1:nmax, function(x) climate_sim_daily_wg[[x]])

    dailyPerformanceCheck(daily.sim = daily_sim_check,
                        daily.obs = daily_obs_check,
                        out.path = out_path,
                        variables = wg_variables[c(1,3,4)],
                        variable.labels = wg_variable_labs[c(1,3,4)],
                        variable.units = wg_variable_units[c(1,3,4)])
  }


  #:::::::::::::::::::::: CLIMATE CHANGE PERTURBATIONS :::::::::::::::::::::::::::

  message(cat("\u2713", "|", "Climate change perturbation matrix created"))


  for (i in 1:length(PARCC)) {

    if(PARCC[[i]]$change_op == "multiply") {

      PARCC[[i]]$mean$begin <- rep(1, 12)
      PARCC[[i]]$var$begin <- rep(1, 12)

    } else {
      PARCC[[i]]$mean$begin <- rep(0, 12)
      PARCC[[i]]$var$begin <- rep(1, 12)
    }

    PARCC[[i]]$mean$steps <-  sapply(1:12, function(m)
      seq(PARCC[[i]]$mean$begin[m], PARCC[[i]]$mean$end[m], length.out = PARCC[[i]]$steps)) %>%
      rbind(rep(NA,12))

    PARCC[[i]]$var$steps <-  sapply(1:12, function(m)
      seq(PARCC[[i]]$var$begin[m], PARCC[[i]]$var$end[m], length.out = PARCC[[i]]$steps)) %>%
      rbind(rep(NA,12))
  }

  # Scenario matrix
  scn_mat <- expand_grid(nvar=1:nmax, par1=PARCC[[1]]$sind, par2=PARCC[[2]]$sind)
  scn_num <- nrow(scn_mat)

  # Loop through each scenario
  perturbs <- list(mean = list(), var = list())
  perturbs_check <- list()

  for (s in 1:scn_num) {

    # Current synthetic realization & parameters
    ns <- scn_mat$nvar[s]
    ps <- c(scn_mat$par1[s], scn_mat$par2[s])

    # New object to store the results
    climate_rlz_new <- climate_sim_daily_wg[[ns]]

    # Loop over all parameters
    for (i in 1:length(PARCC)) {

      # Current perturbation scenario for each variable
      perturb_pars <- list(mean = PARCC[[i]]$mean$steps[ps[i],],
                           var = PARCC[[i]]$var$steps[ps[i],])

      # Loop over all grid cells
      for (x in 1:ngrids) {

        climate_rlz_new[[x]][[PARCC[[i]]$name]] <- quantileMappingDailyTS(
          parms = perturb_pars,
          base.ts = climate_sim_daily_wg[[n]][[x]] %>% dplyr::select(date, value = PARCC[[i]]$name),
          change.operator = PARCC[[i]]$change_op)
      }

    }

    # Write realization to netcdf file
    writeToNetcdf(
      in.path = nc.path,
      in.file = nc.file,
      out.path = paste0(out_path,"netcdf/"),
      out.file = paste0(project_name, "_scn_",sprintf("%02d", s),".nc"),
      out.variables = wg_variables,
      origin.date = sim_date_begin,
      date.length = length(wg_sim_dates),
      dim.names = nc.dimnames,
      gridded.data = climate_sim_daily_wg[[ns]])
  }

  message(cat("\u2713", "|", "Results outputted to netcdf files"))
}










