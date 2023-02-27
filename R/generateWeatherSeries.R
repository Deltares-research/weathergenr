
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
#' @param month.start the first month of the water year (default value is 1).
#' @param evaluate.model logical value indicating wether to save model evaluation plots
#' @param evaluate.grid.num Number of grid cells to be sampled in the evaluation plots
#' @param warm.subset.criteria A list of statistical parameters used for subsetting from the initial annual simulated series
#' @param mc.wet.quantile wet state threshold (quantile value) for markov-chain modeling
#' @param mc.extreme.quantile extremely wet state threshold (quantile value) for markov-chain modeling
#' @param seed a random seed value (nunmeric)
#' @param compute.parallel logical value indicating whether to run (some) functions in parallel
#' @param num.cores Number of cores to be allocated for parallel computing. If left NULL, maximum possible cores minus one is assigned
#' @param dry.spell.change placeholder
#' @param wet.spell.change placeholder
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
  knn.sample.num = 120,
  mc.wet.quantile = 0.3,
  mc.extreme.quantile = 0.8,
  dry.spell.change = rep(1,12),
  wet.spell.change = rep(1,12),
  evaluate.model = FALSE,
  evaluate.grid.num = 20,
  output.path = getwd(),
  seed = sample.int(1e5,1),
  compute.parallel = TRUE,
  num.cores = NULL)

{

  start_time <- Sys.time()

  # Workaround for rlang warning
  wyear <- month <- day <- year <- 0

  if(is.null(variable.labels)) variable.labels <- variable.names
  if(is.null(variable.units)) variable.units <- rep("", length(variable.names))

  # Number of grids
  grids  <- weather.grid$id
  ngrids <- length(grids)

  message(cat(as.character(Sys.time()), "- Random seed: ", seed))

  if(compute.parallel == TRUE) {

    if(is.null(num.cores)) num.cores <- parallel::detectCores()-1
    message(cat(as.character(Sys.time()), "- Parallel mode ( ", num.cores, " cores)"))

  } else {
    message(cat(as.character(Sys.time()), "- Sequential mode."))
  }


  if (!dir.exists(output.path)) {dir.create(output.path)}
  plots_path <- file.path(output.path, "plots")
  if (!dir.exists(plots_path)) {dir.create(plots_path)}

  # PREPARE DATA MATRICES ::::::::::::::::::::::::::::::::::::::::::::::::::::::

  message(cat(as.character(Sys.time()), "- Historical climate data being prepared"))
  message(cat(as.character(Sys.time()), "- Climate variables included:", paste(variable.names, collapse = ', ')))
  message(cat(as.character(Sys.time()), "- Historical data period:", as.character(weather.date[1]), "-", as.character(weather.date[length(weather.date)])))
  message(cat(as.character(Sys.time()), "- Number of climate grids:", ngrids))

  # Historical dates
  year_seq <- as.numeric(format(weather.date,"%Y"))
  year_start <- year_seq[1]
  year_end   <-year_seq[length(year_seq)]

  dates_d <- tibble(year = as.numeric(format(weather.date,"%Y")),
                    wyear = getWaterYear(weather.date, month.start),
                    month = as.numeric(format(weather.date,"%m")),
                    day = as.numeric(format(weather.date,"%d"))) %>%
      filter(wyear >= year_start & wyear <= year_end) %>%
      mutate(date = as.Date(paste(wyear, month, day, sep = "-")), .before=1) %>%
      mutate(dateo = as.Date(paste(year, month, day, sep = "-")), .before=1)

  year.num <- length(unique(dates_d$wyear))
  wyear_index <- which(weather.date %in% dates_d$dateo)

  # Multivariate list of daily climate data
  climate_d <- lapply(1:ngrids, function(i)
    weather.data[[i]][wyear_index,] %>%
    select(all_of(variable.names)) %>%
    mutate(year=dates_d$wyear, .))

  climate_d_aavg <- Reduce(`+`, climate_d) / ngrids

  # Multivariate list of annual climate data
  climate_a <- lapply(1:ngrids, function(i)
        climate_d[[i]] %>% group_by(year) %>%
        summarize(across({{variable.names}}, mean)) %>%
        ungroup() %>% suppressMessages())

  # Area-averaged annual weather series
  climate_a_aavg <- Reduce(`+`, climate_a) / ngrids


  # Simulated dates
  if(is.null(sim.year.num)) sim.year.num <- year.num
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

  #::::::::::: ANNUAL TIME-SERIES GENERATION USING WARM ::::::::::::::::::::::::

  #####  Wavelet analysis on observed annual series
  warm_variable <- climate_a_aavg %>% pull({{warm.variable}})

  # power spectra analysis of historical series
  warm_power <- waveletAnalysis(variable = warm_variable,
      signif.level = warm.signif.level, plot = TRUE, output.path = plots_path)

  # if there is low-frequency signal
  if(length(warm_power$signif_periods) > 0) {

      # wavelet decomposition of historical series
      wavelet_comps <- waveletDecompose(variable = warm_variable,
          signif.periods = warm_power$signif_periods,
          signif.level = warm.signif.level, plot = TRUE, output.path = plots_path)

      message(cat(as.character(Sys.time()), "- Number of significant low-frequency components:", length(wavelet_comps)-1))
      message(cat(as.character(Sys.time()), "- Annual periodicity (yrs):", paste(warm_power$signif_periods, collapse=",")))

      # Simulate annual series of wavelet variable
      message(cat(as.character(Sys.time()), "- Wavelet AR model: simulating", format(warm.sample.num, big.mark=","), " series"))

      sim_annual <- waveletARIMA(wavelet.components = wavelet_comps,
          sim.year.num = sim.year.num, sim.num = warm.sample.num, seed = seed)

    #if there is no low frequency signal
    } else {

      message(cat(as.character(Sys.time()), "- No low-frequency signals detected"))
      message(cat(as.character(Sys.time()), "- Wavelet AR model: simulating", format(warm.sample.num, big.mark=","), "series"))

      # Remove the mean from the component
      MEAN <- mean(warm_variable)

      MODEL <- forecast::auto.arima((warm_variable - MEAN), max.p = 2,max.q = 2,max.P = 0,max.Q = 0,
        stationary = TRUE, seasonal = FALSE)

      INTERCEPT <- ifelse(length(which(names(MODEL$coef)=="intercept")) > 0,
              as.vector(MODEL$coef)[which(names(MODEL$coef)=="intercept")],0)

      sim_annual <- sapply(1:warm.sample.num, function(x) {set.seed(seed+x)
        stats::simulate(MODEL, sim.year.num, sd = sqrt(MODEL$sigma2))}) + INTERCEPT + MEAN

  }

  message(cat(as.character(Sys.time()), "- Filtering stochastic annual series based on default subsetting criteria"))
  # wavelet analysis on simulated series
  sim_power <- sapply(1:warm.sample.num, function(x)
    waveletAnalysis(sim_annual[, x], signif.level = warm.signif.level)$GWS)

  # Define subset range for annual realizations
  if(is.null(warm.subset.criteria)) {

    warm.subset.criteria = list(
        mean = c(0.95,1.05),
        sd = c(0.85,1.15),
        min = c(0.80,1.20),
        max = c(0.80,1.20),
        power = c(0.40,2.60),
        nonsignif.threshold = 0.75)

    if(!length(warm_power$signif_periods > 0)) {
       warm.subset.criteria$power <- NULL
       warm.subset.criteria$nonsignif.threshold <- NULL
    }

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
       output.path = plots_path,
       bounds = warm.subset.criteria,
       seed = seed,
       save.series = FALSE)

  message(cat(as.character(Sys.time()), "-", ncol(sim_annual_sub$subsetted), "stochastic series filtered"))
  message(cat(as.character(Sys.time()), "-", ncol(sim_annual_sub$sampled), "stochastic series sampled"))

  #::::::::::: TEMPORAL & SPATIAL DISSAGGREGATION (knn & mc) :::::::::::::::::::

  resampled_dates <- as_tibble(matrix(0, nrow=nrow(sim_dates_d),
    ncol=realization.num), .name_repair=~paste0("rlz_", 1:realization.num))

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
      TMAX = climate_d_aavg$temp_max,
      TMIN = climate_d_aavg$temp_min,
      START_YEAR_SIM = sim.year.start,
      k1 = n,
      ymax = sim.year.num,
      dates.d = dates_d,
      sim.dates.d = sim_dates_d,
      knn.annual.sample.num = knn.sample.num,
      dry.spell.change = dry.spell.change,
      wet.spell.change = wet.spell.change,
      YEAR_D = year_seq,
      month.start = month.start,
      wet.quantile = mc.wet.quantile,
      extreme.quantile = mc.extreme.quantile,
      seed = seed + n)
  }

  if (compute.parallel) parallel::stopCluster(cl)

  for(x in 1:ncol(resampled_dates)) {
    resampled_dates[,x] <-dates_d$dateo[match(resampled_ini[[x]], dates_d$date)]
  }

  message(cat(as.character(Sys.time()), "- KNN & MCMC modeling completed"))

  utils::write.csv(sim_dates_d$date, file.path(output.path, "sim_dates.csv"), row.names = FALSE)
  utils::write.csv(resampled_dates, file.path(output.path, "resampled_dates.csv"), row.names = FALSE)

  message(cat(as.character(Sys.time()), paste0("- Resampled saved to", output.path)))

  day_order <- sapply(1:realization.num,
    function(n) match(resampled_dates[[n]], dates_d$dateo))

  rlz <- list()
   for (n in 1:realization.num) {
      rlz[[n]] <- lapply(climate_d, function(x)
        x[day_order[,n],] %>% select(-year))
  }


  #::::::::::: MODEL EVALUATION (OPTIONAL) :::::::::::::::::::::::::::::::::::::

    if(evaluate.model) {

    message(cat(as.character(Sys.time()), "- Comparing observed and simulated climate statistics"))

    # Sample evenly from the grid cells
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

    suppressWarnings(
      evaluateWegen(daily.sim = rlz_sample,
                    daily.obs = obs_sample,
                    output.path = plots_path,
                    variables = variable.names,
                    variable.labels = variable.labels,
                    variable.units = variable.units,
                    realization.num = realization.num,
                    wet.quantile = mc.wet.quantile,
                    extreme.quantile = mc.extreme.quantile)
    )

  } else {
     message(cat(as.character(Sys.time()), "- Comparison of historical and simulated climate statistics skipped"))
  }

  message(cat(as.character(Sys.time()), "- Stochastic weather generation completed. See `run.log`"))
  message(cat(as.character(Sys.time()), "- Elapsed time:", Sys.time() - start_time, "mins"))
  #unlink('weathergenr_run.log')


  return(list(resampled = resampled_dates, dates = sim_dates_d$date))

}






