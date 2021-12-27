
#' Simulate gridded weather function
#'
#' Description goes here....
#'
#' @param climate.data list of data frames of daily weather observations per grid cell. Each data frame, columns are weather variables and rows are daily values.
#' @param climate.grid Data frame of grid cells. Each grid cell is assigned an id starting from 1, x and y coordinate index value, and x and y coordinates.
#' @param year.start numeric value indicating the desired starting year of the historical data as input
#' @param year.num numeric value indicating the total number of years of historical data to be considered as input.
#' @param output.path output path for the weather generator results (string)
#' @param variable.names vector of names for the variables to be included in the weather generator
#' @param variable.labels vector of labels for the weather variables (optional). If no values provided, it is labels will be same as the names
#' @param variable.units vector of units for each of the weather variables (optional). If no values provided, a blank vector is used.
#' @param warm.variable the name of the variable for the wavelet autoregressive mode. Default value is precipitation variable.
#' @param warm.signif.level the significance level for the warm model.
#' @param warm.sample.size number of annual sequeces to be generated from the the warm model
#' @param knn.annual.sample.size number of knn years to be sampled
#' @param sim.year.start numeric value indicating the starting year of the generated time-series
#' @param sim.year.num  numeric value indicating the desired total number of years of simulated weather realizations
#' @param realization.num number of natural variability realizations to be generated.
#' @param month.start the first month of the year (default value is 1). Use a value other than 1 for water-year based analyses
#' @param evaluate.model logical value indicating wether to save model evaluation plots
#' @param evaluate.grid.num Number of grid cells to be sampled in the evaluation plots
#' @param apply.delta.changes Logical value indicating weather delta factors to be applied
#' @param delta.precip list of precipitation changes
#' @param delta.temp List of delta changes
#' @param nc.prefix the prefix string to be added to resulting netcdf files
#' @param apply.step.changes logical value to apply transient or step changes
#' @param nc.template template file for the output netcdf
#' @param warm.subset.criteria placeholder
#' @param mc.wet.threshold placeholder
#' @param mc.extreme.quantile placeholder
#'
#' @return
#' @export
#' @import ggplot2
#' @import tibble
#' @import tidyr
#' @import patchwork
#' @import ncdf4
#' @import dplyr
#' @importFrom sf st_as_sf st_sample st_cast st_coordinates
#' @importFrom utils write.csv
#' @importFrom tidyr expand_grid
simulateWeather <- function(
  climate.data = NULL,
  climate.grid = NULL,
  year.start = NULL,
  year.num = NULL,
  month.start = 1,
  variable.names = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  sim.year.start = 2020,
  sim.year.num = 40,
  realization.num = 5,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.size = 20000,
  warm.subset.criteria = NULL,
  knn.annual.sample.size = 100,
  mc.wet.threshold = 0.3,
  mc.extreme.quantile = 0.8,
  evaluate.model = FALSE,
  evaluate.grid.num = 20,
  apply.delta.changes = TRUE,
  apply.step.changes = TRUE,
  delta.precip = NULL,
  delta.temp = NULL,
  output.path = NULL,
  nc.template = NULL,
  nc.prefix = "clim_change_rlz"
  )

 {

  start_time <- Sys.time()

  # Workaround for rlang warning
  wyear <- month <- day <- year <- 0

  stopifnot(is.numeric(year.start))
  stopifnot(is.numeric(year.num))
  stopifnot(is.numeric(sim.year.start))
  stopifnot(is.numeric(sim.year.num))

  # Number of grids
  grids  <- climate.grid$id
  ngrids <- length(grids)

  if(is.null(variable.labels)) variable.labels <- variable.names
  if(is.null(variable.units)) variable.units <- rep("", length(variable.names))

  #browser()
  message(cat("\u2713", "|", "Historical data loaded (spatial resolution:",
    ngrids, "grid cells)"))

  if (!dir.exists(output.path)) {dir.create(output.path)}
  warm_path <- paste0(output.path, "historical/")
  if (!dir.exists(warm_path)) {dir.create(warm_path)}

  hist.date.start <- as.Date(paste0(year.start,"-01-01"))
  date_seq <- hist.date.start-1 + 1:nrow(climate.data[[1]])
  year_seq <- as.numeric(format(date_seq,"%Y"))

  # Historical data matrix tables
  dates_d <- tibble(year = as.numeric(format(date_seq,"%Y")),
          wyear = getWaterYear(date_seq, month.start),
          month = as.numeric(format(date_seq,"%m")),
          day = as.numeric(format(date_seq,"%d"))) %>%
      filter(wyear >= year.start & wyear < year.start + year.num) %>%
      #filter(wyear >= year.start & wyear <= year_seq[length(date_seq)]) %>%
      mutate(date = as.Date(paste(wyear, month, day, sep = "-")), .before=1)

  # Date indices of simulated series
  sim_year_end <- sim.year.start + sim.year.num
  date_sim <- seq(as.Date(paste(sim.year.start,"-1-01",sep="")),
    as.Date(paste(sim_year_end,"-12-31",sep="")), by="day")

  sim_dates_d <- tibble(year = as.numeric(format(date_sim,"%Y")),
      wyear = getWaterYear(date_sim, month.start),
      month = as.numeric(format(date_sim,"%m")),
      day = as.numeric(format(date_sim,"%d"))) %>%
  filter(month!=2 | day!=29) %>%
  filter(wyear >= sim.year.start+1 & wyear <= sim_year_end) %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-")), .before=1)

  year_series <- sim_dates_d$year
  month_series <- sim_dates_d$month
  year_index <- year_series - min(year_series) + 1

  # Prepare tables for daily, monthly, and annual values
  wyear_index <- which(dates_d$date %in% date_seq)
  dates_wy <- tibble(year=dates_d$wyear, month=dates_d$month, day=dates_d$day)

  climate_d <- lapply(1:ngrids, function(i)
    climate.data[[i]][wyear_index,] %>%
    select({{variable.names}}) %>%
    mutate(dates_wy, .))

  climate_a <- lapply(1:ngrids, function(i)
        climate_d[[i]] %>% group_by(year) %>%
        summarize(across({{variable.names}}, mean)) %>%
        ungroup() %>% suppressMessages())

  # Obtain scalar (area-averaged) weather time-series for each variable
  climate_d_aavg <- Reduce(`+`, climate_d) / ngrids
  climate_a_aavg <- Reduce(`+`, climate_a) / ngrids

  #::::::::::: ANNUAL TIME-SERIES GENERATION USING WARM ::::::::::::::::::::::::

  #####  Wavelet analysis on observed annual series
  warm_variable <- climate_a_aavg %>% pull({{warm.variable}})

  # power spectra analysis of historical series
  warm_power <- waveletAnalysis(variable = warm_variable,
    signif.level = warm.signif.level, plot = TRUE, output.path = warm_path)

  # wavelet decomposition of historical series
  wavelet_comps <- waveletDecompose(variable = warm_variable,
        signif.periods = warm_power$signif_periods,
        signif.level = warm.signif.level, plot = TRUE, output.path = warm_path)

  message(cat("\u2713", "|", "Low-frequency components:", length(wavelet_comps)-1,
    paste0("(", sapply(warm_power$signif_periods, function(x) x[1]), " years)")))

  # stochastic simulation of annual series
  sim_annual <- waveletARIMA(wavelet.components = wavelet_comps,
        sim.year.num = sim.year.num, sim.num = warm.sample.size)

  message(cat("\u2713", "|", format(warm.sample.size, big.mark=","),
    "stochastic series simulated with WARM"))

  # wavelet analysis on simulated series
  sim_power <- sapply(1:warm.sample.size, function(x)
    waveletAnalysis(sim_annual[, x], signif.level = warm.signif.level)$GWS)

  # Define subsetting range for annual realizations
  if(is.null(warm.subset.criteria)) {

   warm.subset.criteria = list(
      mean = c(0.95,1.05),
      sd = c(0.90,1.10),
      min = c(0.90,1.15),
      max = c(0.85,1.10),
      power = c(0.50,2.50),
      nonsignif.threshold = 0.75)

    }

  # subsetting from generated warm series
  sim_annual_sub <- waveletARSubset(
       series.obs = warm_variable,
       series.sim = sim_annual,
       power.obs = warm_power$GWS,
       power.sim = sim_power,
       power.period = warm_power$GWS_period,
       power.signif = warm_power$GWS_signif,
       sample.num = realization.num,
       output.path = warm_path,
       bounds = warm.subset.criteria)

  message(cat("\u2713", "|", ncol(sim_annual_sub$subsetted),
    "stochastic series match subsetting criteria"))
  message(cat("\u2713", "|", ncol(sim_annual_sub$sampled),
    "series randomly selected"))

  #::::::::::: TEMPORAL & SPATIAL DISSAGGREGATION (knn & mc) :::::::::::::::::::

  dates_resampled <- lapply(1:realization.num,
    function(n) resampleDates(
        PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
        ANNUAL_PRCP = warm_variable,
        PRCP = climate_d_aavg$precip,
        TEMP = climate_d_aavg$temp,
        START_YEAR_SIM = sim.year.start,
        k1 = n,
        ymax = sim.year.num,
        dates.d = dates_d,
        sim.dates.d = sim_dates_d,
        knn.annual.sample.num = knn.annual.sample.size,
        YEAR_D = year_seq,
        month.start = month.start,
        wet.threshold = mc.wet.threshold,
        extreme.quantile = mc.extreme.quantile
    )
  )

  message(cat("\u2713", "|",
    "Spatially & temporally dissaggregated with knn & mc modeling"))

  dates_resampled_tbl <- bind_cols(dates_resampled,
    .name_repair=~paste0("rlz_", 1:length(dates_resampled)))

  write.csv(dates_resampled_tbl, paste0(warm_path, "dates_resampled.csv"))
  message(cat("\u2713", "|", "Resampled dates saved to",
    paste0(warm_path, "resampled.csv")))

  day_order <- sapply(1:realization.num,
    function(n) match(dates_resampled[[n]], dates_d$date))

  rlz <- list()
   for (n in 1:realization.num) {
      rlz[[n]] <- lapply(climate_d, function(x)
        x[day_order[,n],] %>% select(-year,-month,-day))
  }

  if(evaluate.model) {

    ## Check existing directories and create as needed
    eval_path <- paste0(output.path, "evaluation/")
    if (!dir.exists(eval_path)) dir.create(eval_path)

    sampleGrids <- sf::st_as_sf(climate.grid[,c("x","y")], coords=c("x","y")) %>%
      sf::st_sample(size = min(evaluate.grid.num, ngrids), type="regular") %>%
      sf::st_cast("POINT") %>% sf::st_coordinates() %>% as_tibble() %>%
      left_join(climate.grid[,c("x","y","id")], by = c("X"="x","Y"="y")) %>%
      pull(id)

    rlz_sample <- list()
    for (n in 1:realization.num) {
      rlz_sample[[n]] <- lapply(rlz[[n]][sampleGrids], function(x)
         mutate(x, date = sim_dates_d$date, .before = 1))
    }

    obs_sample <- lapply(climate.data[sampleGrids], function(x)
      dplyr::mutate(x, date = date_seq, .before = 1))

    evaluateWegen(daily.sim = rlz_sample,
       daily.obs = obs_sample,
       output.path = eval_path,
       variables = variable.names,
       variable.labels = variable.labels,
       variable.units = variable.units,
       realization.num = realization.num)

    message(cat("\u2713", "|", "Model evaluation plots saved to:", eval_path))

  }

  #:::::::::::::::::::::: APPLY DELTA CHANGES ::::::::::::::::::::::::::::::::::

  if(!apply.delta.changes) {
      return(dates_resampled_tbl)
  } else {

    # set defaults for the delta factors
    if(is.null(delta.precip$mean$min)) delta.precip$mean$min <- 1
    if(is.null(delta.precip$mean$max)) delta.precip$mean$max <- 1
    if(is.null(delta.precip$var$min)) delta.precip$var$min <- 1
    if(is.null(delta.precip$var$max)) delta.precip$var$max <- 1
    if(is.null(delta.temp$mean$min)) delta.temp$mean$min <- 0
    if(is.null(delta.temp$mean$max)) delta.temp$mean$max <- 0

    # If single value entered, convert to monthly
    if(length(delta.precip$mean$min)==1) rep(delta.precip$mean$min, 12)
    if(length(delta.precip$mean$max)==1) rep(delta.precip$mean$max, 12)
    if(length(delta.precip$var$min)==1) rep(delta.precip$var$min, 12)
    if(length(delta.precip$var$max)==1) rep(delta.precip$var$max, 12)
    if(length(delta.temp$mean$min)==1) rep(delta.temp$mean$min, 12)
    if(length(delta.temp$mean$max)==1) rep(delta.temp$mean$max, 12)

    # Precip and temp change incremental runs
    delta.precip$mean$steps <- sapply(1:12, function(m)
          seq(delta.precip$mean$min[m], delta.precip$mean$max[m],
              length.out = delta.precip$increments))

    delta.precip$var$steps <- sapply(1:12, function(m)
          seq(delta.precip$var$min[m], delta.precip$var$max[m],
              length.out = delta.precip$increments))

    delta.temp$mean$steps <- sapply(1:12, function(m)
          seq(delta.temp$mean$min[m], delta.temp$mean$max[m],
              length.out = delta.temp$increments))

    stresstest_matrix <- tidyr::expand_grid(pind = 1:delta.precip$increments,
          tind = 1:delta.temp$increments) %>%
        mutate(ind = 1:n(), .before = 1)
    smax <- nrow(stresstest_matrix)

     #Create output directory if doesn't exist
    future_path <- paste0(output.path,"future/")
    if (!dir.exists(future_path)) {dir.create(future_path)}

    write.csv(stresstest_matrix, paste0(future_path, "stresstest_matrix.csv"))

    counter <- 0

    # Loop through each scenario
    for (n in 1:realization.num) {

      # Set current weather realization
      rlz_cur <- rlz[[n]]

      # Loop through climate scenarios
      for (s in 1:smax) {

        counter <- counter + 1

        # Current perturbation scenario for each variable
        shift_prcp_mean <- delta.precip$mean$steps[stresstest_matrix$pind[s],]
        shift_prcp_var <- delta.precip$var$steps[stresstest_matrix$pind[s],]
        shift_temp_mean <- delta.temp$mean$steps[stresstest_matrix$tind[s],]

        temp_delta_factors <- sapply(1:12, function(x)
          seq(0, shift_temp_mean[x], length.out = year.num))

        temp_deltas <- sapply(1:length(year_series), function(x)
          temp_delta_factors[year_index[x], month_series[x]])

        # Loop through each grid cell
        for (x in 1:ngrids) {

          # Perturb daily precipitation using quantile mapping
          rlz_cur[[x]]$precip <- quantileMapping(
                  value = rlz_cur[[x]]$precip,
                  mon.ts = month_series,
                  year.ts = year_index,
                  mean.change = shift_prcp_mean,
                  var.change = shift_prcp_var,
                  step.change = apply.step.changes,
                  reltol = 1e-7)

          # Perturb temp, temp_min, and temp_max by delta factors
          rlz_cur[[x]]$temp <- rlz_cur[[x]]$temp + temp_deltas
          rlz_cur[[x]]$temp_min <- rlz_cur[[x]]$temp_min + temp_deltas
          rlz_cur[[x]]$temp_max <- rlz_cur[[x]]$temp_max + temp_deltas

          # Calculate PET from temp, temp_min, temp_max
          rlz_cur[[x]]$pet <- with(rlz_cur[[x]], hargreavesPet(
            months = month_series, temp = temp, tdiff = temp_max - temp_min,
            lat = climate.grid$y[x]))
      }

        # Write to netcdf
        writeNetcdf(
            data = rlz_cur,
            coord.grid = climate.grid,
            output.path = future_path,
            origin.date =  sim_dates_d$date[1],
            calendar.type = "noleap",
            nc.template.file = nc.template,
            nc.compression = 4,
            nc.spatial.ref = "spatial_ref",
            nc.file.prefix = nc.prefix,
            nc.file.suffix = paste0(n,"_", s)
        )




        if (counter < realization.num* smax) {
          cat("\u2059", "|", "Applying climate changes:",
            paste0(counter, "/", realization.num * smax), "\r")
        } else {
          cat("\u2713", "|", "Applying climate changes:",
            paste0(counter, "/", realization.num * smax), "\n")
          message(cat("\u2713", "|", "Elapsed time:",
            (Sys.time() - start_time)/60))
        }

      } # smax close

    } #realization.num close

  } # apply climate change close

} # function close



