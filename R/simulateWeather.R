
#' Simulate gridded weather function
#'
#' @param climate.data placeholder
#' @param climate.grid placeholder
#' @param year.start placeholder
#' @param year.num placeholder
#' @param output.path placeholder
#' @param variable.names placeholder
#' @param variable.labels placeholder
#' @param variable.units placeholder
#' @param warm.variable placeholder
#' @param warm.signif.level placeholder
#' @param warm.sample.size placeholder
#' @param knn.annual.sample.size placeholder
#' @param save.warm.results placeholder
#' @param sim.year.start placeholder
#' @param sim.year.num placeholder
#' @param realization.num placeholder
#' @param month.start placeholder
#' @param evaluate.model placeholder
#' @param evaluate.grid.num placeholder
#' @param apply.climate.changes placeholder
#' @param precip.changes placeholder
#' @param temp.changes placeholder
#' @param output.ncfile.prefix placeholder
#' @param save.scenario.matrix placeholder
#' @param apply.step.changes placeholder
#' @param ... placeholder
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
  output.path = NULL,
  variable.names = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.size = 10000,
  knn.annual.sample.size = 20,
  save.warm.results = TRUE,
  sim.year.start = 2020,
  sim.year.num = 40,
  realization.num = 5,
  month.start = 1,
  evaluate.model = FALSE,
  evaluate.grid.num = 20,
  apply.climate.changes = FALSE,
  precip.changes = NULL,
  temp.changes = NULL,
  save.scenario.matrix = TRUE,
  apply.step.changes = TRUE,
  output.ncfile.prefix = "clim_change_rlz",
  output.ncfile.template = NULL,
  output.ncfile.dimnames = NULL,
  ...)

 {

  start_time <- Sys.time()

  # Number of grids
  grids  <- climate.grid$id
  ngrids <- length(grids)

  if(is.null(variable.labels)) {variable.labels <- variable.names}
  if(is.null(variable.units)) {variable.units <- rep("", length(variable.names))}

  #browser()
  message(cat("\u2713", "|", "Historical data loaded (spatial resolution:", ngrids, "grid cells)"))

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
      filter(wyear >= year.start & wyear <= year_seq[length(date_seq)]) %>%
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
  year_num <- length(unique(year_series))

  # Prepare tables for daily, monthly, and annual values
  wyear_index <- which(dates_d$date %in% date_seq)
  dates_wy <- tibble(year=dates_d$wyear, month=dates_d$month, day=dates_d$day)

  climate_d <- lapply(1:ngrids, function(i)
    climate.data[[i]][wyear_index,] %>%
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
  warm_variable_org <- climate_a_aavg %>% pull({{warm.variable}})
  warm_variable <- warm_variable_org

  # power spectra analysis of historical series
  warm_power <- waveletAnalysis(variable = warm_variable,
    signif.level = warm.signif.level, plot = TRUE, output.path = warm_path)

  # wavelet decomposition of historical series
  wavelet_comps <- waveletDecompose(variable = warm_variable,
        signif.periods = warm_power$signif_periods,
        signif.level = warm.signif.level, plot = TRUE, output.path = warm_path)

  message(cat("\u2713", "|", "Low frequency components:", length(wavelet_comps)-1,
    paste0("(", sapply(warm_power$signif_periods, function(x) x[1]), " years)")))

  # stochastic simulation of annual series
  sim_annual <- waveletARIMA(wavelet.components = wavelet_comps,
        sim.year.num = sim.year.num, sim.num = warm.sample.size)

  message(cat("\u2713", "|", format(warm.sample.size, big.mark=","), "stochastic series simulated with WARM"))

  # wavelet analysis on simulated series
  sim_power <- sapply(1:warm.sample.size, function(x)
    waveletAnalysis(sim_annual[, x], signif.level = warm.signif.level)$GWS)

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
       ...)

  message(cat("\u2713", "|", ncol(sim_annual_sub$subsetted), "stochastic series match subsetting criteria"))
  message(cat("\u2713", "|", ncol(sim_annual_sub$sampled), "series randomly selected"))

  # Sample size in KNN_ANNUAL sampling
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
        month.start = month.start)
  )

  message(cat("\u2713", "|", "Spatially & temporally dissaggregated with knn & mc modeling"))

  dates_resampled_tbl <- bind_cols(dates_resampled,
    .name_repair=~paste0("rlz_", 1:length(dates_resampled)))

  write.csv(dates_resampled_tbl, paste0(warm_path, "dates_resampled.csv"))
  message(cat("\u2713", "|", "Resampled dates outputted to csv file"))

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

    sampleGrids <- sf::st_as_sf(climate.grid[,c("x","y")], coords = c("x","y")) %>%
      sf::st_sample(size = min(evaluate.grid.num, ngrids), type = "regular") %>%
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

  #:::::::::::::::::::::: CLIMATE CHANGE SETTINGS ::::::::::::::::::::::::::::::

  if(!apply.climate.changes) {
      return(dates_resampled_tbl)
  } else {

    precip.changes$mean$steps <- sapply(1:12, function(m)
          seq(precip.changes$mean$min[m], precip.changes$mean$max[m],
              length.out = precip.changes$increments))

    precip.changes$var$steps <- sapply(1:12, function(m)
          seq(precip.changes$var$min[m], precip.changes$var$max[m],
              length.out = precip.changes$increments))

    temp.changes$mean$steps <- sapply(1:12, function(m)
          seq(temp.changes$mean$min[m], temp.changes$mean$max[m],
              length.out = temp.changes$increments))

    scn_mat_index <- tidyr::expand_grid(precip_ind = 1:precip.changes$increments,
          temp_ind = 1:temp.changes$increments) %>%
        mutate(ind = 1:n(), .before = 1)
    smax <- nrow(scn_mat_index)

     #Create output directory if doesn't exist
    future_path <- paste0(output.path,"future/")
    if (!dir.exists(future_path)) {dir.create(future_path)}

    if(isTRUE(save.scenario.matrix)) {
      write.csv(x = scn_mat_index, file = paste0(future_path, "scenario_matrix.csv"))
    }

    counter <- 0

    # Loop through each scenario
    for (n in 1:realization.num) {

      # Set current weather realization
      rlz_cur <- rlz[[n]]

      # Loop through climate scenarios
      for (s in 1:smax) {

        counter <- counter + 1

        # Current perturbation scenario for each variable
        perturb_precip_mean <- precip.changes$mean$steps[scn_mat_index$precip_ind[s],]
        perturb_precip_var <- precip.changes$var$steps[scn_mat_index$precip_ind[s],]
        perturb_temp_mean <- temp.changes$mean$steps[scn_mat_index$temp_ind[s],]

        temp_delta_factors <- sapply(1:12, function(x)
          seq(0, perturb_temp_mean[x], length.out = year_num))

        temp_deltas <- sapply(1:length(year_series), function(x)
          temp_delta_factors[year_index[x], month_series[x]])

        # Loop through each grid cell
        for (x in 1:ngrids) {

          # Perturb daily precipitation using quantile mapping
          rlz_cur[[x]]$precip <- quantileMapping(value = rlz_cur[[x]]$precip,
                  mon.ts = month_series, year.ts = year_index,
                  mean.change = perturb_precip_mean,
                  var.change = perturb_precip_var,
                  step.change = apply.step.changes,
                  reltol = 1e-5)

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
            nc.temp = output.ncfile.template,
            data = rlz_cur,
            coord.grid = climate.grid,
            output.path = future_path,
            nc.dimnames = output.ncfile.dimnames,
            origin.date =  sim_dates_d$date[1],
            calendar.type = "no leap",
            variables = c(variable.names, "pet"),
            variable.units = c(variable.units, "mm/day"),
            file.prefix = output.ncfile.prefix,
            file.suffix = paste0(n,"_", s)
        )

        if (counter < realization.num* smax) {
          cat("\u2059", "|", "Applying climate changes:", paste0(counter, "/", realization.num * smax), "\r")
        } else {
          cat("\u2713", "|", "Applying climate changes:", paste0(counter, "/", realization.num * smax), "\n")
          message(cat("\u2713", "|", "Elapsed time:", (Sys.time() - start_time)/60))
        }

      } # smax close


    } #realization.num close

  } # apply climate change close

} # function close



