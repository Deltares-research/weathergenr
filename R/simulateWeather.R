
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
  rmax = 5000,
  nmax = 5,
  nc.dimnames = NULL,
  validate = FALSE,
  ...)

 {

  # Number of grids
  grids  <- climate_tidy$id
  ngrids <- length(grids)

  #browser()
  message(cat("\u2713", "|", "Historical data loaded. Data has", ngrids, "grids."))

  if (!dir.exists(output.dir)) {dir.create(output.dir)}
  out_path_warm <- paste0(output.dir, "warm/")
  if (!dir.exists(out_path_warm)) {dir.create(out_path_warm)}

  ##### Set historical gridded data to be inputted to wg (remove leap days!)
  noleap_index <- with(climate_tidy$data[[1]], which(!(month(date) == 2 & day(date) == 29)))
  climate_obs_daily <- lapply(climate_tidy$data, function(x) x[noleap_index, ])

  #Reference data date variables
  wg_dates <- climate_obs_daily[[1]]$date
  # Month order for water year
  wy_months <- setdiff(c(month(wg_dates[1]):12, 1:(month(wg_dates[1])-1)),0)
  # Date template table with water-year adjustment
  wg_dates_adjusted <- tibble(year = year(wg_dates) - min(year(wg_dates)) + 1,
      mon = month(wg_dates), day = day(wg_dates)) %>%
      mutate(mon = as.numeric(as.character(factor(mon, levels = wy_months, labels= 1:12))))

  WATER_YEAR_A = unique(year(wg_dates))
  month_list = wy_months
  water_year_start = WATER_YEAR_A[1]
  water_year_end = WATER_YEAR_A[length(WATER_YEAR_A)]

  # Daily format
  DATE_D = wg_dates
  MONTH_D = month(wg_dates)
  DAY_D = day(wg_dates)
  WATER_YEAR_D = year(wg_dates)
  YEAR_D = WATER_YEAR_D
  MONTH_DAY_D = as.matrix(data.frame(MONTH_D=MONTH_D, DAY_D=DAY_D))

  # Date indices of simulated series
  START_YEAR_SIM <- 2020
  END_YEAR_SIM <- START_YEAR_SIM + ymax
  DATE_SIM <- seq(as.Date(paste(START_YEAR_SIM,"-1-01",sep="")),as.Date(paste(END_YEAR_SIM,"-12-31",sep="")),by="day")
  DAY_SIM <- as.numeric(format(DATE_SIM,"%d"))
  MONTH_SIM <- as.numeric(format(DATE_SIM,"%m"))
  YEAR_SIM <- as.numeric(format(DATE_SIM,"%Y"))
  no_leap <- which(MONTH_SIM!=2 | DAY_SIM!= 29)
  DAY_SIM <- DAY_SIM[no_leap]
  MONTH_SIM <- MONTH_SIM[no_leap]
  YEAR_SIM <- YEAR_SIM[no_leap]

  WATER_YEAR_SIM <- YEAR_SIM
  #if (water_yr_change) {WATER_YEAR_SIM[which(MONTH_SIM>=month_list[1])] <- WATER_YEAR_SIM[which(MONTH_SIM>=month_list[1])] + 1}
  WATER_YEAR_LOCATIONS_SIM <- which(WATER_YEAR_SIM>=(START_YEAR_SIM+1) & WATER_YEAR_SIM<=END_YEAR_SIM)
  WATER_YEAR_SIM <- WATER_YEAR_SIM[WATER_YEAR_LOCATIONS_SIM]
  YEAR_SIM <- YEAR_SIM[WATER_YEAR_LOCATIONS_SIM]
  MONTH_SIM <- MONTH_SIM[WATER_YEAR_LOCATIONS_SIM]
  DAY_SIM <- DAY_SIM[WATER_YEAR_LOCATIONS_SIM]
  DATE_SIM <- as.Date(paste(WATER_YEAR_SIM,MONTH_SIM,DAY_SIM,sep="-"))

  SIM_LENGTH <- length(DATE_SIM)
  DATE_M_SIM <- subset(DATE_SIM,DAY_SIM==1)
  YEAR_M_SIM <- as.numeric(format(DATE_M_SIM,"%Y"))
  MONTH_M_SIM <- as.numeric(format(DATE_M_SIM,"%m"))
  DATE_A_SIM <- subset(DATE_M_SIM,MONTH_M_SIM==month_list[1])
  WATER_YEAR_A_SIM <- as.numeric(format(DATE_A_SIM,"%Y"))

  sim_dates <- tibble(date = seq.Date(wg_dates[1], wg_dates[1]-1 + years(ymax), by = "day")) %>%
    filter(!(month(date)==2 & day(date)==29)) %>% pull(date)

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
  warm_power <- waveletAnalysis(variable = warm_annual, variable.unit = "mm",
    signif.level = ssig, plot = TRUE, out.path = out_path_warm)

  ##### WAVELET DECOMPOSITION OF HISTORICAL SERIES
  wavelet_comps <- waveletDecompose(variable = warm_annual,
        signif.periods = warm_power$signif_periods,
        signif.level = ssig, plot = TRUE, out.path = out_path_warm)

  message(cat("\u2713", "|", "Wavelet analysis:", length(wavelet_comps)-1, "low frequency components defined"))

  ##### Stochastic simuluation of historical annual series
  sim_annual <- waveletAR(wavelet.comps = wavelet_comps,
        num.years = ymax, num.realizations = rmax)

  message(cat("\u2713", "|", rmax, "stochastic annual series generated"))

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
       out.path = out_path_warm,
       ...)
       #mean.bounds = c(0.95,1.05),
       #sdev.bounds = c(0.95,1.10),
       #max.bounds  = c(0.95,1.10),
       #min.bounds  = c(0.90,1.05),
       #power.bounds = c(0.8,2.5),
       #nonsig.threshold = 0.9)

  # Subsetted realizations of annual simulated time-series
  sim_annual_final <- sim_annual[,c(1,2,3,4,5)]
  sim_annual_final <- sim_annual_sub$sampled

  # Plot simulated warm-series
  df1 <- sim_annual_sub$sampled %>%
    as_tibble(.name_repair = ~paste0("rlz",1:nmax)) %>%
    mutate(x = 1:ymax) %>%
    gather(key = variable, value=y, -x) %>%
    mutate(y = y * 365)

  df2 <- tibble(x=1:length(warm_annual_org), y = warm_annual_org*365)

  (ggplot(df1, aes(x = x, y = y)) +
    theme_light() +
    geom_line(aes(y = y, group = variable, color = variable), alpha = 0.6) +
    geom_line(aes(y=y), data = df2, color = "black", size = 1) +
    scale_x_continuous(limits = c(0,ymax), breaks = seq(0,ymax, 5)) +
    #scale_y_continuous(limits = c(1500, 3500), breaks = seq(1500,3500,500)) +
    guides(color = "none") +
    scale_color_brewer(palette = "PuOr") +
    labs(y = "Precip (mm)", x = "Year index")) %>%
  ggsave(filename = paste0(out_path_warm, "warm_simulated_series.png"),.,
    height = 5, width = 10)

  message(cat("\u2713", "|", ncol(sim_annual_sub$sampled), "annual traces selected"))

  #::::MCMC MODELING COMES HERE! ::::::::::::;;;;;;;;;;::::::::::::::::::::::::::::::

  # Sample size in KNN_ANNUAL sampling
  kk <- max(round(sqrt(length(warm_annual)),0), round(length(warm_annual),0)*.5)

  rlz_daily_index <- lapply(1:nmax, function(n)
    DAILY_WEATHER_GENERATOR(k1 = n, ymax = ymax,
      PRCP_FINAL_ANNUAL_SIM = sim_annual_final[, n], ANNUAL_PRCP = warm_annual,
      WATER_YEAR_A = WATER_YEAR_A, WATER_YEAR_D = WATER_YEAR_D,
      PRCP = climate_daily_wg_aavg$precip, TEMP = climate_daily_wg_aavg$temp, DATE_D = DATE_D, MONTH_D = MONTH_D, YEAR_D = YEAR_D,
      MONTH_DAY_D = MONTH_DAY_D, month_list = month_list, water_year_start = water_year_start,
      water_year_end = water_year_end, y_sample_size = 20, SIM_LENGTH = length(sim_dates),
      kk = kk, MONTH_SIM = MONTH_SIM, WATER_YEAR_SIM = WATER_YEAR_SIM, START_YEAR_SIM = START_YEAR_SIM,
      DAY_SIM = DAY_SIM)
  )

  sim_daily_wg <- list()

  for (n in 1:nmax) {

    day_order <- match(rlz_daily_index[[n]], wg_dates)
    sim_daily_wg[[n]] <- lapply(climate_daily_wg, function(x)
      x[day_order,] %>% mutate(date = sim_dates, .before = 1))
  }

  message(cat("\u2713", "|", "Spatial & temporal dissaggregation through markov-chain knn completed."))

  #:::::::::::::::::::::: PERFORMANCE EVALUATION OF DAILY SIMULATED DATA :::::::

  if(isTRUE(validate)) {

    ## Check existing directories and create as needed

    out_path_performance <- paste0(output.dir, "performance/")
    if (!dir.exists(out_path_performance)) {dir.create(out_path_performance)}

    sample_ngrid <- min(20, ngrids)

    sampleGrids <- sf::st_as_sf(climate_tidy[,c("x","y")], coords = c("x","y")) %>%
      sf::st_sample(size = sample_ngrid, type = "regular") %>%
      sf::st_cast("POINT") %>% sf::st_coordinates() %>%
      as_tibble() %>%
      left_join(climate_tidy[,c("x","y","id")], by = c("X"="x","Y"="y")) %>%
      pull(id)

    sim_daily_wg_sample <- lapply(1:nmax, function(x) sim_daily_wg[[x]][sampleGrids])

    dailyPerformance(daily.sim = sim_daily_wg_sample,
                        daily.obs = climate_obs_daily[sampleGrids],
                        out.path = out_path_performance,
                        variables = wg.vars,
                        variable.labels = wg.var.labs,
                        variable.units = wg.var.units,
                        nmax = nmax)

    message(cat("\u2713", "|", "Validation plots saved to output folder"))
  }

  #::::::::::::::::::::: OUTPUT HISTORICAL REALIZATIONS TO FILE ::::::::::::::::


  out_path_hist <- paste0(output.dir, "historical/")
  if (!dir.exists(out_path_hist)) {dir.create(out_path_hist)}

  # Dimensions in the outputted variable (order matters!)
  dim_ord <- names(nc_data$nc_dimensions)

  ncout_dims <- list()
  ncout_dims[[nc.dimnames$time]] <- ncdf4::ncdim_def(name = nc.dimnames$time,
      units = paste0("days since ", format(round(as.POSIXct(wg.date.begin), units = "day"),
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
    ncout_file <- ncdf4::nc_create(paste0(out_path_hist, proj.name,"_hist_rlz_",n, ".nc"),
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
   message(cat("\u2713", "|", nmax, "Daily gridded weather realizations saved to output folder"))

}
