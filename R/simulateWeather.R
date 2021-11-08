
#' Generate daily, multivariate griddded weather series
#'
#' @param output.path path to save the results
#' @param hist.climate list of data frames for historical daily weather data with variable
#' @param hist.date.start Starting year for the historical data to be used in the analysis, e.g., 1980
#' @param hist.date.end Ending year for the historical data to be used in the analysis, e.g., 2000
#' @param variables Names of weather variables to be included, e.g., precip, temp, temp_min, temp_max.
#' @param variable.labels Long names of weather variables. If not provided, same as variable names.
#' @param variable.units Measurment units of the weather variables
#' @param warm.variable Variable to be used for the wavelet analysis
#' @param warm.signif.level Statistical significance level used in the wavelet analysis.
#' @param warm.sample.num Placeholder.
#' @param sim.year.start Placeholder.
#' @param sim.year.num Placeholder.
#' @param realization.num Placeholder.
#' @param check.stats Placeholder.
#' @param check.grid.num Placeholder.

#' @param month.start Placeholder.
#' @param ... Placeholder
#' @param grid.coords Placeholder
#' @param knn.annual.sample.num Placeholder
#' @param return.date.indices Placeholder
#' @param save.to.netcdf Placeholder
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
simulateWeather <- function(
  output.path = NULL,
  hist.climate = NULL,
  grid.coords = NULL,
  hist.year.start = NULL,
  hist.year.num = NULL,
  variables = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.num = 20000,
  knn.annual.sample.num = 20,
  sim.year.num = 40,
  sim.year.start = 2020,
  realization.num = 5,
  check.stats = FALSE,
  month.start = 5,
  check.grid.num = 50,
  return.date.indices = TRUE,
  save.to.netcdf = FALSE,
  ...)

 {

  # Number of grids
  ngrids <- length(hist.climate)

  if(is.null(variable.labels)) {variable.labels <- variables}
  if(is.null(variable.units)) {variable.units <- rep("", length(variables))}

  #browser()
  message(cat("\u2713", "|", "Historical data loaded. Data has", ngrids, "grids"))

  if (!dir.exists(output.path)) {dir.create(output.path)}
  warm_path <- paste0(output.path, "warm/")
  if (!dir.exists(warm_path)) {dir.create(warm_path)}

  hist.date.start <- as.Date(paste0(hist.year.start,"-01-01"))
  hist.date.end <- as.Date(paste0(hist.year.start+hist.year.num,"-01-01"))

  date_seq <- hist.date.start-1 + 1:nrow(hist.climate[[1]])
  year_seq <- as.numeric(format(date_seq,"%Y"))

  # Historical data matrix tables
  dates_d <- tibble(year = as.numeric(format(date_seq,"%Y")),
          wyear = getWaterYear(date_seq, month.start),
          month = as.numeric(format(date_seq,"%m")),
          day = as.numeric(format(date_seq,"%d"))) %>%
      filter(wyear >= hist.year.start & wyear <= hist.year.start+hist.year.num) %>%
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

  # Prepare tables for daily, monthly, and annual values
  wyear_index <- which(dates_d$date %in% date_seq)
  dates_wy <- tibble(year=dates_d$wyear, month=dates_d$month, day=dates_d$day)
  climate_d <- lapply(1:ngrids, function(i)
    hist.climate[[i]][wyear_index,] %>%
    bind_cols(dates_wy, .))

  climate_a <- lapply(1:ngrids, function(i)
        climate_d[[i]] %>% group_by(year) %>%
        summarize(across({{variables}}, mean)) %>%
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

  message(cat("\u2713", "|", "Wavelet analysis:", length(wavelet_comps)-1,
    "low frequency components defined"))

  # stochastic simulation of annual series
  sim_annual <- waveletARIMA(wavelet.components = wavelet_comps,
        sim.year.num = sim.year.num, sim.num = warm.sample.num)

  message(cat("\u2713", "|", warm.sample.num, "stochastic annual series generated"))

  # wavelet analysis on simulated series
  sim_power <- sapply(1:warm.sample.num, function(x)
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

  # Subsetted realizations of annual simulated time-series
  sim_annual_final <- sim_annual_sub$sampled

  message(cat("\u2713", "|", ncol(sim_annual_final), "annual traces selected"))

  # Sample size in KNN_ANNUAL sampling
  dates_resampled <- lapply(1:realization.num,
    function(n) resampleDates(
        PRCP_FINAL_ANNUAL_SIM = sim_annual_final[, n],
        ANNUAL_PRCP = warm_variable,
        PRCP = climate_d_aavg$precip,
        TEMP = climate_d_aavg$temp,
        START_YEAR_SIM = sim.year.start,
        k1 = n,
        ymax = sim.year.num,
        dates.d = dates_d,
        sim.dates.d = sim_dates_d,
        knn.annual.sample.num = knn.annual.sample.num,
        YEAR_D = year_seq,
        month.start = month.start)
  )



  message(cat("\u2713", "|", "Spatial & temporal dissaggregation completed"))

   day_order <- sapply(1:realization.num,
    function(n) match(dates_resampled[[n]], date_seq))

  rlz <- list()
   for (n in 1:realization.num) {
      rlz[[n]] <- lapply(climate_d, function(x)
        x[day_order[,n],] %>% select(-year,-month,-day))
  }


  if(check.stats) {

    ## Check existing directories and create as needed
    if (!dir.exists(paste0(output.path, "check/"))) {
              dir.create(paste0(output.path, "check/"))}

    sampleGrids <- sf::st_as_sf(grid.coords[,c("x","y")], coords = c("x","y")) %>%
      sf::st_sample(size = min(check.grid.num, ngrids), type = "regular") %>%
      sf::st_cast("POINT") %>% sf::st_coordinates() %>%
      as_tibble() %>%
      left_join(grid.coords[,c("x","y","id")], by = c("X"="x","Y"="y")) %>%
      pull(id)

    rlz_sample <- list()
    for (n in 1:realization.num) {
      rlz_sample[[n]] <- lapply(rlz[[n]][sampleGrids], function(x)
         mutate(x, date = sim_dates_d$date, .before = 1))
    }

    obs_sample <- lapply(hist.climate[sampleGrids], function(x)
      dplyr::mutate(x, date = date_seq, .before = 1))

    wegenCheck(daily.sim = rlz_sample,
       daily.obs = obs_sample,
       output.path = paste0(output.path, "check/"),
       variables = variables,
       variable.labels = variable.labels,
       variable.units = variable.units,
       realization.num = realization.num)

    message(cat("\u2713", "|", "control plots saved"))

  }

  if(save.to.netcdf) {

    for (n in 1:realization.num) {

      writeNetcdf(
        data = rlz[[n]],
        coord.grid = nc_data$coords,
        output.path = paste0(out_path,"historical/"),
        nc.dimnames = list(x = "lon", y = "lat", time = "time"),
        origin.date = as.Date(paste0(sim.year.start,"-01-01")),
        calendar.type = "no leap",
        variables = variables,
        variable.units = variable.units,
        file.suffix = n)
    }

    message(cat("\u2713", "|", "Save results to netcdf files"))

  }

  if(return.date.indices) {
    message(cat("\u2713", "|", "output date orders"))
    return(day_order)

  } else {
    message(cat("\u2713", "|", "output multivariate series"))
    return(rlz)
  }

}



