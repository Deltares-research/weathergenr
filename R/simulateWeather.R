
#' Main wrapper function for daily weather generation
#'
#' @param output.path Placeholder.
#' @param hist.climate Placeholder.
#' @param hist.date.start Placeholder.
#' @param variables Placeholder.
#' @param variable.labels Placeholder.
#' @param variable.units Placeholder.
#' @param warm.variable Placeholder.
#' @param warm.signif.level Placeholder.
#' @param ymax Placeholder.
#' @param rmax Placeholder.
#' @param nmax Placeholder.
#' @param validate Placeholder.
#' @param validation.grid.num Placeholder.
#' @param sim.year.start Placeholder.
#' @param month.list Placeholder.
#' @param ... Placeholder.
#'
#' @return
#' @export
#' @import ggplot2 tibble tidyr patchwork ncdf4 dplyr
#' @importFrom sf st_as_sf st_sample st_cast st_coordinates
simulateWeather <- function(
  output.path = NULL,
  hist.climate = NULL,
  hist.date.start = NULL,
  variables = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  warm.variable = NULL,
  warm.signif.level = 0.90,
  knn.annual.sample.size = 20,
  ymax = 40,
  rmax = 5000,
  nmax = 5,
  validate = TRUE,
  sim.year.start = 2020,
  month.list = 1:12,
  validation.grid.num = 50,
  return.sampled.date.indices = TRUE,
  ...)

 {

  # Number of grids
  grids  <- hist.climate$id
  ngrids <- length(grids)

  #browser()
  message(cat("\u2713", "|", "Historical data loaded. Data has", ngrids, "grids."))

  if (!dir.exists(output.path)) {dir.create(output.path)}
  warm_path <- paste0(output.path, "warm/")
  if (!dir.exists(warm_path)) {dir.create(warm_path)}

  date_seq <- hist.date.start-1 + 1:nrow(hist.climate$data[[1]])
  year_seq <- as.numeric(format(date_seq,"%Y"))

  water_year_start <- year_seq[1]
  water_year_end <- year_seq[length(year_seq)]

  # Historical data matrix tables
  dates_d <- tibble(year = as.numeric(format(date_seq,"%Y")),
          wyear = getWaterYear(date=date_seq, month.list=month.list),
          month = as.numeric(format(date_seq,"%m")),
          day = as.numeric(format(date_seq,"%d"))) %>%
      filter(wyear >= water_year_start & wyear <= water_year_end) %>%
      mutate(date = as.Date(paste(wyear, month, day, sep = "-")), .before=1)
  dates_m <- dates_d %>% dplyr::filter(day == 1)

  dates_a <- dates_m %>% dplyr::filter(month == month.list[1])
  wyear_index <- which(dates_d$date %in% date_seq)

  # Historical data for the correct period
  climate_obs_ini <- lapply(hist.climate$data, function(x) x)
  climate_obs_daily <- lapply(climate_obs_ini, function(x) x[wyear_index, ])

  # Date indices of simulated series
  sim.end.year <- sim.year.start + ymax

  date_sim <- seq(as.Date(paste(sim.year.start,"-1-01",sep="")),
    as.Date(paste(sim.end.year,"-12-31",sep="")), by="day")

  sim_dates_d <- tibble(year = as.numeric(format(date_sim,"%Y")),
          wyear = getWaterYear(date=date_sim, month.list=month.list),
          month = as.numeric(format(date_sim,"%m")),
          day = as.numeric(format(date_sim,"%d"))) %>%
      filter(month!=2 | day!=29) %>%
      filter(wyear >= (sim.year.start+1) & wyear <= sim.end.year) %>%
      mutate(date = as.Date(paste(wyear, month, day, sep = "-")), .before=1)

  sim_dates_m <- sim_dates_d %>% dplyr::filter(day == 1)
  sim_dates_a <- sim_dates_m %>% dplyr::filter(month == month.list[1])
  sim_wyear_index <- which(sim_dates_d$date %in% date_sim)

  # Prepare tables for daily, monthly, and annual values
  dates_wy <- tibble(year=dates_d$wyear, month=dates_d$month, day=dates_d$day)

  climate_d <- lapply(1:ngrids, function(i)
        bind_cols(dates_wy, climate_obs_daily[[i]]) %>%
        arrange(year, month, day))

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

  ############
  # Normality tests come here.........
  ###########

  warm_variable <- warm_variable_org

  ####  Power spectra of observed annual series
  warm_power <- waveletAnalysis(variable = warm_variable, variable.unit = "mm",
    signif.level = warm.signif.level, plot = TRUE, out.path = warm_path)

  ##### WAVELET DECOMPOSITION OF HISTORICAL SERIES
  wavelet_comps <- waveletDecompose(variable = warm_variable,
        signif.periods = warm_power$signif_periods,
        signif.level = warm.signif.level, plot = TRUE, out.path = warm_path)

  message(cat("\u2713", "|", "Wavelet analysis:", length(wavelet_comps)-1,
    "low frequency components defined"))

  ##### Stochastic simuluation of historical annual series
  sim_annual <- waveletAR(wavelet.comps = wavelet_comps,
        num.years = ymax, num.realizations = rmax)

  message(cat("\u2713", "|", rmax, "stochastic annual series generated"))

  # Wavelet analysis on simulated series
  sim_power <- sapply(1:rmax, function(x)
    waveletAnalysis(sim_annual[, x], signif.level = warm.signif.level)$GWS)

  ### Choose which parameters to consider for filtering
  sim_annual_sub <- waveletARSubset(
       series.obs = warm_variable,
       series.sim = sim_annual,
       power.obs = warm_power$GWS,
       power.sim = sim_power,
       power.period = warm_power$GWS_period,
       power.signif = warm_power$GWS_signif,
       nmax = nmax,
       out.path = warm_path,
       ...)

  # Subsetted realizations of annual simulated time-series
  sim_annual_final <- sim_annual_sub$sampled

  # Plot simulated warm-series
  df1 <- sim_annual_final %>%
    as_tibble(.name_repair = ~paste0("rlz",1:nmax)) %>%
    mutate(x = 1:ymax) %>%
    gather(key = variable, value=y, -x) %>%
    mutate(y = y * 365)

  df2 <- tibble(x=1:length(warm_variable_org), y = warm_variable_org*365)

  (ggplot(df1, aes(x = x, y = y)) +
    theme_light(base_size = 12) +
    geom_line(aes(y = y, group = variable, color = variable), alpha = 0.6) +
    geom_line(aes(y=y), data = df2, color = "black", size = 1) +
    scale_x_continuous(limits = c(0,ymax), breaks = seq(0,ymax, 5)) +
    #scale_y_continuous(limits = c(1500, 3500), breaks = seq(1500,3500,500)) +
    guides(color = "none") +
    scale_color_brewer(palette = "PuOr") +
    labs(y = "Precip (mm)", x = "Year index")) %>%
  ggsave(filename = paste0(warm_path, "warm_simulated_series.png"),.,
    height = 5, width = 10)

  message(cat("\u2713", "|", ncol(sim_annual_final), "annual traces selected"))

  #::::MCMC MODELING COMES HERE! ::::::::::::;;;;;;;;;;::::::::::::::::::::::::::::::

  # Sample size in KNN_ANNUAL sampling
  kk <- round(max(round(sqrt(length(warm_variable)),0), round(length(warm_variable),0)*.5))

  dates_resampled <- lapply(1:nmax, function(n)
    resampleDates(
      ANNUAL_PRCP = warm_variable,
      PRCP_FINAL_ANNUAL_SIM = sim_annual_final[, n],
      WATER_YEAR_A = dates_a$wyear,
      WATER_YEAR_D = dates_d$wyear,
      PRCP = climate_d_aavg$precip,
      TEMP = climate_d_aavg$temp,
      DATE_D = dates_d$date,
      MONTH_D = dates_d$month,
      YEAR_D = year_seq,
      MONTH_DAY_D = dates_d[,c("month","day")],
      month_list = month.list,
      water_year_start = water_year_start,
      water_year_end = water_year_end,
      knn.annual.sample.size = knn.annual.sample.size,
      k1 = n,
      ymax = ymax,
      SIM_LENGTH = length(sim_dates_d$date),
      kk = kk,
      MONTH_SIM = sim_dates_d$month,
      WATER_YEAR_SIM = sim_dates_d$wyear,
      START_YEAR_SIM = sim.year.start,
      DAY_SIM = sim_dates_d$day)
  )

  message(cat("\u2713", "|", "Spatial & temporal dissaggregation completed."))

  sim_daily_wg <- list()
  for (n in 1:nmax) {

      day_order <- match(dates_resampled[[n]], date_seq)
      sim_daily_wg[[n]] <- lapply(climate_d, function(x)
        x[day_order,] %>% mutate(date = sim_dates_d$date, .before = 1) %>%
          select(-year,-month,-day))
  }


  if(validate) {

    ## Check existing directories and create as needed
    out_path_performance <- paste0(output.path, "performance/")
    if (!dir.exists(out_path_performance)) {dir.create(out_path_performance)}

    sampleGrids <- sf::st_as_sf(hist.climate[,c("x","y")], coords = c("x","y")) %>%
      sf::st_sample(size = min(validation.grid.num, ngrids), type = "regular") %>%
      sf::st_cast("POINT") %>% sf::st_coordinates() %>%
      as_tibble() %>%
      left_join(hist.climate[,c("x","y","id")], by = c("X"="x","Y"="y")) %>%
      pull(id)

    sim_daily_sample <- lapply(1:nmax, function(x) sim_daily_wg[[x]][sampleGrids])
    hist_daily_sample <- lapply(climate_obs_ini[sampleGrids], function(x)
      dplyr::mutate(x, date = date_seq, .before = 1))

    wgPerformance(daily.sim = sim_daily_sample,
       daily.obs = hist_daily_sample,
       out.path = out_path_performance,
       variables = variables,
       variable.labels = variable.labels,
       variable.units = variable.units,
       nmax = nmax)

    message(cat("\u2713", "|", "Validation plots saved to output folder"))

  }

  if(save.to.netcdf) {

    message(cat("\u2713", "|", "saving to netcdf files"))

    for (n in 1:nmax) {

      writeNetcdf(
        data = sim_daily_wg[[n]],
        coord.grid = nc_data$tidy_data,
        output.path = paste0(out_path,"historical/"),
        nc.dimensions = nc_data$nc_dimensions,
        nc.dimnames = nc_dimnames,
        origin.date = sim_origin_date,
        calendar.type = "no leap",
        variables = wg_variables,
        variable.units = wg_variable_units,
        file.suffix = n)
    }
  }

  if(return.sampled.date.indices) {
    message(cat("\u2713", "|", "output date orders"))
    return(day_order)

  } else {
    message(cat("\u2713", "|", "output multivariate series"))
    return(sim_daily_wg)
  }

}



