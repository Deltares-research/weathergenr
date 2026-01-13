#' @title Generate Synthetic Gridded Daily Weather Series
#'
#' @description
#' Generate stochastic gridded daily weather by coupling an annual-scale
#' low-frequency generator with daily-scale resampling and persistence logic.
#' The workflow combines:
#' \itemize{
#'   \item Wavelet Autoregressive Modeling (WARM) on an annual aggregate of
#'         \code{warm.variable} to simulate low-frequency variability,
#'   \item annual K-nearest-neighbor (KNN) matching to select historical analogue years,
#'   \item a three-state (dry, wet, extreme) daily Markov chain to control spell persistence,
#'   \item daily KNN resampling of precipitation and temperature anomalies.
#' }
#'
#' The simulation regime is inferred from \code{month.start}:
#' \itemize{
#'   \item \code{month.start == 1}: calendar-year simulation,
#'   \item \code{month.start != 1}: water-year simulation starting at \code{month.start}.
#' }
#'
#' @details
#' \strong{Calendar handling (robust Gregorian input):}
#' Inputs may be on the Gregorian calendar and may include Feb 29. Internally,
#' the function enforces a 365-day calendar by removing Feb 29 from \code{weather.date}
#' and dropping the corresponding rows from every element of \code{weather.data}.
#' All downstream processing and outputs therefore use a 365-day calendar.
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item \strong{Historical preprocessing:} enforce a 365-day calendar, compute
#'         calendar/water-year indices, and build a historical dates table used for
#'         KNN matching and Markov chain sequencing.
#'   \item \strong{Annual WARM simulation:} aggregate historical daily data to annual
#'         means (by water year if applicable), run wavelet analysis on the annual
#'         series, simulate candidate annual traces via \code{\link{wavelet_arima}},
#'         and subset traces using \code{\link{filter_warm_simulations}}.
#'   \item \strong{Daily disaggregation:} for each realization, call
#'         \code{resample_weather_dates()} to generate a simulated sequence of
#'         historical analogue dates (in the internal 365-day calendar).
#'   \item \strong{Output construction:} map resampled internal dates back to
#'         historical observation dates (\code{dateo}) and return simulated dates.
#' }
#'
#' @param weather.data Named list of data frames (one per grid cell) with observed
#'   daily weather. Each data frame must have one row per day corresponding to
#'   \code{weather.date}, and columns for all \code{variables}.
#' @param weather.grid Data frame describing grid cells. Must contain columns
#'   \code{id}, \code{xind}, \code{yind}, \code{x}, \code{y}.
#' @param weather.date Vector of class \code{Date} corresponding to rows of each
#'   \code{weather.data[[i]]}. May include Feb 29; Feb 29 is removed internally
#'   along with the corresponding rows in \code{weather.data}.
#' @param variables Character vector of daily variables to simulate
#'   (e.g., \code{c("precip","temp")}).
#' @param sim.year.num Integer number of years to simulate. If \code{NULL},
#'   defaults to the number of unique historical simulation-years (calendar year
#'   or water year depending on \code{month.start}).
#' @param sim.year.start Integer first simulation year (calendar year if
#'   \code{month.start == 1}, otherwise first water year).
#' @param month.start Integer in \code{1:12}. First month of the simulation year.
#' @param realization.num Integer number of realizations to generate.
#' @param warm.variable Character name of the annual driver variable used in WARM
#'   (must be in \code{variables}).
#' @param warm.signif.level Numeric in (0,1). Wavelet significance level for
#'   retaining low-frequency components in WARM.
#' @param warm.sample.num Integer number of candidate annual traces to generate
#'   before subsetting.
#' @param knn.sample.num Integer number of historical years sampled in annual KNN.
#' @param mc.wet.quantile Numeric in (0,1). Wet-day threshold quantile.
#' @param mc.extreme.quantile Numeric in (0,1). Extreme-day threshold quantile.
#' @param dry.spell.change Numeric length-12 vector of monthly dry-spell adjustment factors.
#' @param wet.spell.change Numeric length-12 vector of monthly wet-spell adjustment factors.
#' @param output.path Character directory for diagnostics/CSV outputs.
#' @param seed Optional integer seed for reproducibility.
#' @param compute.parallel Logical; if \code{TRUE}, run disaggregation in parallel.
#' @param num.cores Optional integer number of cores for parallel execution.
#'
#' @return A list with:
#' \describe{
#'   \item{resampled}{Tibble/data.frame with one column per realization. Each
#'     column contains the historical \code{dateo} values selected as analogues
#'     for each simulated day.}
#'   \item{dates}{Vector of \code{Date} giving the simulated daily time axis
#'     (Feb 29 excluded).}
#' }
#'
#' @export
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @import dplyr
#' @import logger
#' @importFrom dplyr mutate
#' @title Generate Synthetic Gridded Daily Weather Series
#'
#' @description
#' Generate stochastic gridded daily weather by coupling an annual-scale
#' low-frequency generator with daily-scale resampling and persistence logic.
#' The workflow combines:
#' \itemize{
#'   \item Wavelet Autoregressive Modeling (WARM) on an annual aggregate of
#'         \code{warm.variable} to simulate low-frequency variability,
#'   \item annual K-nearest-neighbor (KNN) matching to select historical analogue years,
#'   \item a three-state (dry, wet, extreme) daily Markov chain to control spell persistence,
#'   \item daily KNN resampling of precipitation and temperature anomalies.
#' }
#'
#' The simulation regime is inferred from \code{month.start}:
#' \itemize{
#'   \item \code{month.start == 1}: calendar-year simulation,
#'   \item \code{month.start != 1}: water-year simulation starting at \code{month.start}.
#' }
#'
#' @details
#' \strong{Calendar handling (robust Gregorian input):}
#' Inputs may be on the Gregorian calendar and may include Feb 29. Internally,
#' the function enforces a 365-day calendar by removing Feb 29 from \code{weather.date}
#' and dropping the corresponding rows from every element of \code{weather.data}.
#' All downstream processing and outputs therefore use a 365-day calendar.
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item \strong{Historical preprocessing:} enforce a 365-day calendar, compute
#'         calendar/water-year indices, and build a historical dates table used for
#'         KNN matching and Markov chain sequencing.
#'   \item \strong{Annual WARM simulation:} aggregate historical daily data to annual
#'         means (by water year if applicable), run wavelet analysis on the annual
#'         series, simulate candidate annual traces via \code{\link{wavelet_arima}},
#'         and subset traces using \code{\link{filter_warm_simulations}}.
#'   \item \strong{Daily disaggregation:} for each realization, call
#'         \code{resample_weather_dates()} to generate a simulated sequence of
#'         historical analogue dates (in the internal 365-day calendar).
#'   \item \strong{Output construction:} map resampled internal dates back to
#'         historical observation dates (\code{dateo}) and return simulated dates.
#' }
#'
#' @param weather.data Named list of data frames (one per grid cell) with observed
#'   daily weather. Each data frame must have one row per day corresponding to
#'   \code{weather.date}, and columns for all \code{variables}.
#' @param weather.grid Data frame describing grid cells. Must contain columns
#'   \code{xind}, \code{yind}, \code{x}, \code{y}. If \code{id} is missing, it is
#'   created as a sequential index.
#' @param weather.date Vector of class \code{Date} corresponding to rows of each
#'   \code{weather.data[[i]]}. May include Feb 29; Feb 29 is removed internally
#'   along with the corresponding rows in \code{weather.data}.
#' @param variables Character vector of daily variables to simulate
#'   (e.g., \code{c("precip","temp")}).
#' @param sim.year.num Integer number of years to simulate. If \code{NULL},
#'   defaults to the number of unique historical simulation-years (calendar year
#'   or water year depending on \code{month.start}).
#' @param sim.year.start Integer first simulation year (calendar year if
#'   \code{month.start == 1}, otherwise first water year).
#' @param month.start Integer in \code{1:12}. First month of the simulation year.
#' @param realization.num Integer number of realizations to generate.
#' @param warm.variable Character name of the annual driver variable used in WARM
#'   (must be in \code{variables}).
#' @param warm.signif.level Numeric in (0,1). Wavelet significance level for
#'   retaining low-frequency components in WARM.
#' @param warm.sample.num Integer number of candidate annual traces to generate
#'   before subsetting.
#' @param knn.sample.num Integer number of historical years sampled in annual KNN.
#' @param mc.wet.quantile Numeric in (0,1). Wet-day threshold quantile.
#' @param mc.extreme.quantile Numeric in (0,1). Extreme-day threshold quantile.
#' @param dry.spell.change Numeric length-12 vector of monthly dry-spell adjustment factors.
#' @param wet.spell.change Numeric length-12 vector of monthly wet-spell adjustment factors.
#' @param output.path Character directory for diagnostics/CSV outputs.
#' @param seed Optional integer seed for reproducibility.
#' @param compute.parallel Logical; if \code{TRUE}, run disaggregation in parallel.
#' @param num.cores Optional integer number of cores for parallel execution.
#' @param verbose Logical; if \code{TRUE}, emits progress logs via \code{.log_info()}.
#'
#' @return A list with:
#' \describe{
#'   \item{resampled}{Tibble/data.frame with one column per realization. Each
#'     column contains the historical \code{dateo} values selected as analogues
#'     for each simulated day.}
#'   \item{dates}{Vector of \code{Date} giving the simulated daily time axis
#'     (Feb 29 excluded).}
#' }
#'
#' @export
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @import dplyr
#' @import logger
#' @importFrom dplyr mutate
generate_weather_series <- function(
    weather.data = NULL,
    weather.grid = NULL,
    weather.date = NULL,
    variables = NULL,
    sim.year.num = NULL,
    sim.year.start = 2020,
    month.start = 1,
    realization.num = 5,
    warm.variable = "precip",
    warm.signif.level = 0.90,
    warm.sample.num = 5000,
    knn.sample.num = 120,
    mc.wet.quantile = 0.3,
    mc.extreme.quantile = 0.8,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    output.path = tempdir(),
    seed = NULL,
    compute.parallel = FALSE,
    num.cores = NULL,
    verbose = FALSE
) {
  start_time <- Sys.time()
  verbose <- isTRUE(verbose)

  .log <- function(msg, tag = NULL) {
    .log_info(msg, verbose = verbose, tag = tag)
  }

  # ---------------------------------------------------------------------------
  # Setup
  # ---------------------------------------------------------------------------
  dir.create(output.path, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  stopifnot(
    "weather.data must be a list" = is.list(weather.data),
    "weather.grid must be a data frame" = is.data.frame(weather.grid),
    "weather.date must be a Date vector" = inherits(weather.date, "Date"),
    "weather.data must be non-empty" = length(weather.data) >= 1L,
    "weather.date length must match data rows" = length(weather.date) == nrow(weather.data[[1]]),
    "month.start must be 1-12" = month.start %in% 1:12,
    "variables must be non-empty" = length(variables) > 0,
    "variables must exist in weather.data" = all(variables %in% names(weather.data[[1]])),
    "dry.spell.change must have length 12" = length(dry.spell.change) == 12,
    "wet.spell.change must have length 12" = length(wet.spell.change) == 12,
    "mc.wet.quantile must be in (0,1)" = is.numeric(mc.wet.quantile) && mc.wet.quantile > 0 && mc.wet.quantile < 1,
    "mc.extreme.quantile must be in (0,1)" = is.numeric(mc.extreme.quantile) && mc.extreme.quantile > 0 && mc.extreme.quantile < 1,
    "mc.extreme.quantile > mc.wet.quantile" = mc.extreme.quantile > mc.wet.quantile,
    "warm.variable must be in variables" = warm.variable %in% variables,
    "weather.grid must have required columns" = all(c("id", "xind", "yind", "x", "y") %in% names(weather.grid))
  )

  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE/FALSE.", call. = FALSE)
  }

  if (length(unique(vapply(weather.data, nrow, integer(1)))) != 1L) {
    stop("All grid cells must have the same number of observations", call. = FALSE)
  }

  if (!("id" %in% names(weather.grid))) {
    weather.grid$id <- seq_len(nrow(weather.grid))
  }

  # ---------------------------------------------------------------------------
  # RNG management
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  warm_seed <- sample.int(.Machine$integer.max, 1L)
  daily_seed <- sample.int(.Machine$integer.max, 1L)

  # ---------------------------------------------------------------------------
  # Parallel computing setup
  # ---------------------------------------------------------------------------
  if (isTRUE(compute.parallel)) {
    if (is.null(num.cores)) num.cores <- max(1L, parallel::detectCores() - 1L)

    cl <- parallel::makeCluster(num.cores)
    doParallel::registerDoParallel(cl)

    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = daily_seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    .log("Starting in parallel mode", tag = "INIT")
    .log("Number of cores: {num.cores}", tag = "INIT")
  } else {
    .log("Starting in sequential mode", tag = "INIT")
  }

  # ---------------------------------------------------------------------------
  # Initialization logging
  # ---------------------------------------------------------------------------
  n_grids <- length(weather.data)
  .log("Randomization seed: {seed}", tag = "INIT")
  .log("Climate variables: {paste(variables, collapse = ', ')}", tag = "INIT")
  .log("Total number of grids: {n_grids}", tag = "INIT")
  .log("Historical period: {weather.date[1]} to {weather.date[length(weather.date)]}", tag = "INIT")

  # ---------------------------------------------------------------------------
  # Calendar normalization (enforce 365-day calendar)
  # ---------------------------------------------------------------------------
  leap_idx <- find_leap_days(weather.date)
  if (!is.null(leap_idx) && length(leap_idx) > 0L) {
    weather.date_old <- weather.date
    weather.date <- weather.date[-leap_idx]

    weather.data <- lapply(weather.data, function(df) {
      if (!is.data.frame(df)) stop("Each element of weather.data must be a data.frame.", call. = FALSE)
      if (nrow(df) != length(weather.date_old)) {
        stop("Each weather.data[[i]] must have nrow equal to length(weather.date).", call. = FALSE)
      }
      df[-leap_idx, , drop = FALSE]
    })

    .log("Dropped {length(leap_idx)} row(s): 365-day calendar enforced", tag = "INIT")
  }

  # ---------------------------------------------------------------------------
  # Historical date table (internal time axis)
  # ---------------------------------------------------------------------------
  his_date <- weather.date
  his_wyear <- get_water_year(his_date, month.start)

  dates_d <- tibble::tibble(
    dateo = his_date,
    year  = as.integer(format(his_date, "%Y")),
    wyear = his_wyear,
    month = as.integer(format(his_date, "%m")),
    day   = as.integer(format(his_date, "%d"))
  )

  dates_d <- dates_d |>
    dplyr::mutate(
      date = if (month.start == 1) {
        dateo
      } else {
        as.Date(sprintf("%04d-%02d-%02d", wyear, month, day))
      }
    )

  if (anyNA(dates_d$date)) stop("Internal error: dates_d$date contains NA.", call. = FALSE)
  if (anyDuplicated(dates_d$date)) {
    stop("Internal error: dates_d$date contains duplicates; mapping is ambiguous.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Restrict to complete years (365 days) and keep longest contiguous block
  # ---------------------------------------------------------------------------
  wyear_counts <- as.data.frame(table(dates_d$wyear), stringsAsFactors = FALSE)
  names(wyear_counts) <- c("wyear", "n_days")
  wyear_counts$wyear <- as.integer(wyear_counts$wyear)

  full_wyears <- wyear_counts$wyear[wyear_counts$n_days == 365L]
  if (length(full_wyears) == 0L) {
    stop(
      "No complete simulation-years found in historical record (need 365 days per year after Feb 29 removal).",
      call. = FALSE
    )
  }

  full_wyears <- sort(full_wyears)
  runs <- split(full_wyears, cumsum(c(1L, diff(full_wyears) != 1L)))
  full_wyears <- runs[[which.max(lengths(runs))]]

  dates_d <- dplyr::filter(dates_d, wyear %in% full_wyears)

  wyear_idx <- match(dates_d$dateo, weather.date)
  if (anyNA(wyear_idx)) stop("Internal error: dates_d$dateo did not match weather.date.", call. = FALSE)

  .log("Using {length(full_wyears)} complete year(s) ({min(full_wyears)}-{max(full_wyears)})", tag = "INIT")

  # ---------------------------------------------------------------------------
  # Area-averaged daily and annual climate series
  # ---------------------------------------------------------------------------
  if (n_grids == 1L) {
    climate_d_aavg <- weather.data[[1]][wyear_idx, variables, drop = FALSE]
  } else {
    n_days <- length(wyear_idx)
    n_vars <- length(variables)

    climate_d_aavg_mat <- matrix(0, nrow = n_days, ncol = n_vars)
    colnames(climate_d_aavg_mat) <- variables

    for (i in seq_len(n_grids)) {
      grid_data <- as.matrix(weather.data[[i]][wyear_idx, variables, drop = FALSE])
      climate_d_aavg_mat <- climate_d_aavg_mat + grid_data
    }

    climate_d_aavg <- as.data.frame(climate_d_aavg_mat / n_grids)
  }

  climate_d_aavg$wyear <- dates_d$wyear

  climate_a_aavg <- climate_d_aavg |>
    dplyr::group_by(wyear) |>
    dplyr::summarize(dplyr::across(dplyr::all_of(variables), mean), .groups = "drop")

  # ---------------------------------------------------------------------------
  # Simulation date table (internal 365-day calendar)
  # ---------------------------------------------------------------------------
  if (is.null(sim.year.num)) sim.year.num <- length(unique(dates_d$wyear))
  sim.year.num <- as.integer(sim.year.num)

  sim_year_end <- sim.year.start + sim.year.num
  sim_date_start <- as.Date(sprintf("%04d-%02d-01", sim.year.start, month.start))
  sim_date_end <- as.Date(sprintf("%04d-%02d-01", sim_year_end, month.start)) - 1

  sim_date_ini <- seq.Date(sim_date_start, sim_date_end, by = "day")
  sim_date_ini <- sim_date_ini[format(sim_date_ini, "%m-%d") != "02-29"]

  sim_dates_d <- tibble::tibble(
    dateo = sim_date_ini,
    year  = as.integer(format(sim_date_ini, "%Y")),
    wyear = get_water_year(sim_date_ini, month.start),
    month = as.integer(format(sim_date_ini, "%m")),
    day   = as.integer(format(sim_date_ini, "%d"))
  ) |>
    dplyr::mutate(
      date = if (month.start == 1) {
        dateo
      } else {
        as.Date(sprintf("%04d-%02d-%02d", wyear, month, day))
      }
    )

  if (anyNA(sim_dates_d$date)) stop("Internal error: sim_dates_d$date contains NA.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # Annual time-series generation (WARM)
  # ---------------------------------------------------------------------------
  .log("Running annual WARM simulation", tag = "WARM")

  detrend <- TRUE
  warm_variable <- climate_a_aavg %>% dplyr::pull({{ warm.variable }})

  warm_power <- wavelet_spectral_analysis(
    variable = warm_variable,
    signif.level = warm.signif.level,
    period.lower.limit = 2,
    detrend = detrend,
    mode = "complete"
  )

  p <- plot_wavelet_spectra(
    variable = warm_variable,
    variable.year = climate_a_aavg$wyear,
    period = warm_power$gws_period,
    POWER = warm_power$power,
    GWS = warm_power$gws_unmasked,
    GWS_signif = warm_power$gws_signif_unmasked,
    coi = warm_power$coi,
    sigm = warm_power$sigm
  )

  tryCatch(
    ggplot2::ggsave(file.path(output.path, "global_wavelet_power_spectrum.png"), p, width = 8, height = 5),
    error = function(e) logger::log_warn("Failed to save wavelet spectrum plot: {e$message}")
  )

  if (any(!is.na(warm_power$signif_periods))) {
    .log("Significant low-frequency components detected: {length(warm_power$comps_names) - 1}", tag = "WARM")
    .log(
      "Detected periodicity (years): {paste(round(warm_power$gws_period[warm_power$signif_periods], 2), collapse = ', ')}",
      tag = "WARM"
    )
  } else {
    .log("No low-frequency signals detected", tag = "WARM")
  }

  sim_annual <- wavelet_arima(
    wavelet.components = warm_power$comps,
    sim.year.num = sim.year.num,
    sim.num = warm.sample.num,
    seed = warm_seed,
    match.variance = TRUE
  )

  sim_annual_sub <- filter_warm_simulations(
    series.obs = warm_variable,
    series.sim = sim_annual,
    sample.num = realization.num,
    seed = warm_seed + 1L,
    padding = TRUE,
    bounds = list(),
    relax.order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
    make.plots = TRUE,
    wavelet.pars = list(
      signif.level = warm.signif.level,
      noise.type = "red",
      period.lower.limit = 2,
      detrend = detrend
    ),
    verbose = verbose
  )

  tryCatch(
    {
      ggplot2::ggsave(file.path(output.path, "warm_annual_series.png"), sim_annual_sub$plots[[1]], width = 8, height = 5)
      ggplot2::ggsave(file.path(output.path, "warm_annual_statistics.png"), sim_annual_sub$plots[[2]], width = 8, height = 5)
      ggplot2::ggsave(file.path(output.path, "warm_annual_wavelet.png"), sim_annual_sub$plots[[3]], width = 8, height = 5)
    },
    error = function(e) logger::log_warn("Failed to save warm plots: {e$message}")
  )

  # ---------------------------------------------------------------------------
  # Daily disaggregation (KNN + Markov chain) -> resampled historical dates
  # ---------------------------------------------------------------------------
  .log("Starting daily weather simulation using KNN + Markov Chain scheme", tag = "KNN")

  if (isTRUE(compute.parallel)) {
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }

  resampled_ini <- foreach::foreach(n = seq_len(realization.num)) %d% {
    resample_weather_dates(
      annual_prcp_sim = sim_annual_sub$sampled[, n],
      annual_prcp_obs = warm_variable,
      daily_prcp_obs = climate_d_aavg$precip,
      daily_temp_obs = climate_d_aavg$temp,
      sim_start_year = sim.year.start,
      realization_id = n,
      n_sim_years = sim.year.num,
      obs_dates_table = dates_d,
      sim_dates_table = sim_dates_d,
      annual_knn_n = knn.sample.num,
      year_start_month = month.start,
      wet_q = mc.wet.quantile,
      extreme_q = mc.extreme.quantile,
      dry_spell_factor = dry.spell.change,
      wet_spell_factor = wet.spell.change,
      seed = daily_seed + n
    )
  }

  # ---------------------------------------------------------------------------
  # Map internal simulated dates -> original historical observation dates (dateo)
  # ---------------------------------------------------------------------------
  resampled_dates <- vector("list", realization.num)
  names(resampled_dates) <- paste0("rlz_", seq_len(realization.num))

  for (i in seq_len(realization.num)) {
    resampled_dates[[i]] <- dates_d$dateo[match(resampled_ini[[i]], dates_d$date)]
  }
  resampled_dates <- tibble::as_tibble(resampled_dates)

  # ---------------------------------------------------------------------------
  # Outputs
  # ---------------------------------------------------------------------------
  tryCatch(
    utils::write.csv(sim_dates_d$date, file.path(output.path, "sim_dates.csv"), row.names = FALSE),
    error = function(e) logger::log_error("Failed to write sim_dates.csv: {e$message}")
  )

  tryCatch(
    utils::write.csv(resampled_dates, file.path(output.path, "resampled_dates.csv"), row.names = FALSE),
    error = function(e) logger::log_error("Failed to write resampled_dates.csv: {e$message}")
  )

  .log("Results saved to: `{output.path}`", tag = "DONE")
  .log("Simulation complete. Elapsed time: {Sys.time() - start_time} secs", tag = "DONE")

  list(resampled = resampled_dates, dates = sim_dates_d$dateo)
}
