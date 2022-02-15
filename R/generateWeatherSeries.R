
#' Simulate gridded weather function
#'
#' Description goes here....
#'
#' @param weather.data list of data frames of daily weather observations per grid cell. Each data frame, columns are weather variables and rows are daily values.
#' @param weather.date a vector of dates matching the weather.data
#' @param weather.grid Data frame of grid cells. Each grid cell is assigned an id starting from 1, x and y coordinate index value, and x and y coordinates.
#' @param output.path output path for the weather generator results (string)
#' @param variable.names vector of names for the variables to be included in the weather generator
#' @param variable.labels vector of labels for the weather variables (optional). If no values provided, it is labels will be same as the names
#' @param variable.units vector of units for each of the weather variables (optional). If no values provided, a blank vector is used.
#' @param warm.variable the name of the variable for the wavelet autoregressive mode. Default value is precipitation variable.
#' @param warm.signif.level the significance level for the warm model.
#' @param warm.sample.num number of annual sequeces to be generated from the the warm model
#' @param knn.sample.num number of knn years to be sampled
#' @param sim.year.start numeric value indicating the starting year of the generated time-series
#' @param sim.year.num  numeric value indicating the desired total number of years of simulated weather realizations
#' @param realization.num number of natural variability realizations to be generated.
#' @param month.start the first month of the year (default value is 1). Use a value other than 1 for water-year based analyses
#' @param evaluate.model logical value indicating wether to save model evaluation plots
#' @param evaluate.grid.num Number of grid cells to be sampled in the evaluation plots
#' @param warm.subset.criteria A list of statistical parameters used for subsetting from the initial annual simulated series
#' @param mc.wet.quantile wet state threshold (quantile value) for markov-chain modeling
#' @param mc.extreme.quantile extremely wet state threshold (quantile value) for markov-chain modeling
#' @param seed a random seed value (nunmeric)
#' @param compute.parallel logical value indicating whether to run (some) functions in parallel
#' @param num.cores Number of cores to be allocated for parallel computing. If left NULL, maximum possible cores minus one is assigned
#'
#'
#' @return
#' @export
#' @import ggplot2
#' @import tibble
#' @import tidyr
#' @import patchwork
#' @import dplyr
generateWeatherSeries <- function(
  weather.data = NULL,
  weather.grid = NULL,
  weather.date = NULL,
  variable.names = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  sim.year.num = NULL,
  sim.year.start = 2020,
  month.start = 1,
  realization.num = 5,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.num = 5000,
  warm.subset.criteria = NULL,
  knn.sample.num = 100,
  mc.wet.quantile = 0.3,
  mc.extreme.quantile = 0.8,
  evaluate.model = FALSE,
  evaluate.grid.num = 20,
  output.path = getwd(),
  seed = NULL,
  compute.parallel = TRUE,
  num.cores = NULL)

 {

  start_time <- Sys.time()

  # Workaround for rlang warning
  wyear <- month <- day <- year <- 0

  if(is.null(seed)) seed <- sample.int(1e5,1)
  if(is.null(variable.labels)) variable.labels <- variable.names
  if(is.null(variable.units)) variable.units <- rep("", length(variable.names))

  # Number of grids
  grids  <- weather.grid$id
  ngrids <- length(grids)

  if(compute.parallel == TRUE) {

    if(is.null(num.cores)) num.cores <- parallel::detectCores()-1
    message(cat("\u2713", "|", paste0("Stochastic time-series generation in parallel mode (", num.cores, " cores)")))

  } else {
    message(cat("\u2713", "|", "Stochastic time-series generation in sequential mode"))
  }

  #browser()
  message(cat("\u2713", "|",
      paste0("Input weather data: ", ngrids, " grid cells, ",
    length(variable.names), " variables with a length of ", length(weather.date), " days")
  ))

  if (!dir.exists(output.path)) {dir.create(output.path)}
  warm_path <- paste0(output.path, "historical/")
  if (!dir.exists(warm_path)) {dir.create(warm_path)}

  # PREPARE DATA MATRICES ::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # Historical dates
  year_seq <- as.numeric(format(weather.date,"%Y"))
  year_start <- year_seq[1]
  year_end   <-year_seq[length(year_seq)]

  dates_d <- tibble(year = as.numeric(format(weather.date,"%Y")),
          wyear = getWaterYear(weather.date, month.start),
          month = as.numeric(format(weather.date,"%m")),
          day = as.numeric(format(weather.date,"%d"))) %>%
      filter(wyear >= year_start & wyear <= year_end) %>%
      mutate(date = as.Date(paste(wyear, month, day, sep = "-")), .before=1)

  year.num <- length(unique(dates_d$wyear))
  if(is.null(sim.year.num)) sim.year.num <- year.num

  wyear_index <- which(dates_d$date %in% weather.date)
  dates_wy <- tibble(year=dates_d$wyear, month=dates_d$month,
    day=dates_d$day)

  # Simulated dates
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

  # Daily and annual climate data
  climate_d <- lapply(1:ngrids, function(i)
    weather.data[[i]][wyear_index,] %>%
    select(all_of(variable.names)) %>%
    mutate(dates_wy, .))
  climate_d_aavg <- Reduce(`+`, climate_d) / ngrids

  climate_a <- lapply(1:ngrids, function(i)
        climate_d[[i]] %>% group_by(year) %>%
        summarize(across({{variable.names}}, mean)) %>%
        ungroup() %>% suppressMessages())
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

  message(cat("\u2713", "|", "Number of low-frequency components:", length(wavelet_comps)-1,
    paste0("(periodicity: ", sapply(warm_power$signif_periods, function(x) x[1]), " years)")))

  # Simulate annual series of wavelet variable
  sim_annual <- waveletARIMA(wavelet.components = wavelet_comps,
        sim.year.num = sim.year.num, sim.num = warm.sample.num, seed = seed)

  message(cat("\u2713", "|", format(warm.sample.num, big.mark=","),
    "stochastic series simulated with WARM"))

  # wavelet analysis on simulated series
  sim_power <- sapply(1:warm.sample.num, function(x)
    waveletAnalysis(sim_annual[, x], signif.level = warm.signif.level)$GWS)

  # Define subsetting range for annual realizations
  if(is.null(warm.subset.criteria)) {

     warm.subset.criteria = list(
        mean = c(0.95,1.05),
        sd = c(0.85,1.15),
        min = c(0.80,1.20),
        max = c(0.80,1.20),
        power = c(0.40,2.60),
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
       bounds = warm.subset.criteria,
       seed = seed)

  message(cat("\u2713", "|", ncol(sim_annual_sub$subsetted),
    "stochastic series match subsetting criteria"))

  message(cat("\u2713", "|", ncol(sim_annual_sub$sampled),
    "series sampled"))


  #::::::::::: TEMPORAL & SPATIAL DISSAGGREGATION (knn & mc) :::::::::::::::::::

  resampled_dates <- as_tibble(matrix(0, nrow=nrow(sim_dates_d), ncol=realization.num),
    .name_repair=~paste0("rlz_", 1:realization.num))

  if (compute.parallel) {

    cl <- parallel::makeCluster(num.cores)
    doParallel::registerDoParallel(cl)

    `%d%` <- foreach::`%dopar%`

  } else {

    `%d%` <- foreach::`%do%`

  }

  resampled_ini <- foreach::foreach(n=seq_len(realization.num)) %d% {

    weathergenr::resampleDates(
      PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
      ANNUAL_PRCP = warm_variable,
      PRCP = climate_d_aavg$precip,
      TEMP = climate_d_aavg$temp,
      START_YEAR_SIM = sim.year.start,
      k1 = n,
      ymax = sim.year.num,
      dates.d = dates_d,
      sim.dates.d = sim_dates_d,
      knn.annual.sample.num = knn.sample.num,
      YEAR_D = year_seq,
      month.start = month.start,
      wet.quantile = mc.wet.quantile,
      extreme.quantile = mc.extreme.quantile,
      seed = seed + n)
  }

  if (compute.parallel) parallel::stopCluster(cl)

  for(x in 1:ncol(resampled_dates)) {resampled_dates[,x] <- resampled_ini[[x]]}

  message(cat("\u2713", "|",
    "Spatial/temporal dissaggregation with knn & markov chain modeling"))

  utils::write.csv(sim_dates_d$date, paste0(warm_path, "sim_dates.csv"), row.names = FALSE)
  utils::write.csv(resampled_dates, paste0(warm_path, "resampled_dates.csv"), row.names = FALSE)

  message(cat("\u2713", "|", "Resampled dates saved as resampled_dates.csv"))
  day_order <- sapply(1:realization.num,
    function(n) match(resampled_dates[[n]], dates_d$date))

  rlz <- list()
   for (n in 1:realization.num) {
      rlz[[n]] <- lapply(climate_d, function(x)
        x[day_order[,n],] %>% select(-year,-month,-day))
  }


  #::::::::::: MODEL EVALUATION (OPTIONAL) :::::::::::::::::::::::::::::::::::::

  if(evaluate.model) {

    ## Check existing directories and create as needed
    eval_path <- paste0(output.path, "evaluation/")
    if (!dir.exists(eval_path)) dir.create(eval_path)

    sampleGrids <- sf::st_as_sf(weather.grid[,c("x","y")], coords=c("x","y")) %>%
      sf::st_sample(size = min(evaluate.grid.num, ngrids), type="regular") %>%
      sf::st_cast("POINT") %>% sf::st_coordinates() %>% as_tibble() %>%
      left_join(weather.grid[,c("x","y","id")], by = c("X"="x","Y"="y")) %>%
      pull(id)

    rlz_sample <- list()
    for (n in 1:realization.num) {
      rlz_sample[[n]] <- lapply(rlz[[n]][sampleGrids], function(x)
         mutate(x, date = sim_dates_d$date, .before = 1))
    }

    obs_sample <- lapply(weather.data[sampleGrids], function(x)
      dplyr::mutate(x, date = weather.date, .before = 1))

    suppressWarnings(evaluateWegen(daily.sim = rlz_sample,
       daily.obs = obs_sample,
       output.path = eval_path,
       variables = variable.names,
       variable.labels = variable.labels,
       variable.units = variable.units,
       realization.num = realization.num)
    )


    message(cat("\u2713", "|", "Model evaluation plots saved to:", eval_path))

  }

  return(list(resampled = resampled_dates, dates = sim_dates_d$date))

}
