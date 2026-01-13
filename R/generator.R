#' Generate Synthetic Gridded Daily Weather Series
#'
#' @description
#' Generates stochastic gridded daily weather by coupling an annual-scale
#' low-frequency generator with daily-scale resampling and persistence logic.
#' The workflow combines:
#' \itemize{
#'   \item Wavelet Autoregressive Modeling (WARM) on an annual aggregate of
#'     \code{warm_var} to simulate low-frequency variability,
#'   \item annual K-nearest-neighbor (KNN) matching to select historical analogue years,
#'   \item a three-state (dry, wet, extreme) daily Markov chain to control spell persistence,
#'   \item daily KNN resampling of precipitation and temperature anomalies.
#' }
#'
#' The simulation-year definition is inferred from \code{year_start_month}:
#' \itemize{
#'   \item \code{year_start_month == 1}: calendar-year simulation,
#'   \item \code{year_start_month != 1}: water-year simulation starting in \code{year_start_month}.
#' }
#'
#' @details
#' \strong{Calendar handling (robust Gregorian input):}
#' Inputs may be on the Gregorian calendar and may include Feb 29. Internally,
#' the function enforces a 365-day calendar by removing Feb 29 from \code{obs_dates}
#' and dropping the corresponding rows from every element of \code{obs_data}.
#' All downstream processing and outputs therefore use a 365-day calendar.
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item \strong{Historical preprocessing:} Enforce a 365-day calendar, compute
#'     calendar/water-year indices, and build a historical dates table used for
#'     KNN matching and Markov chain sequencing.
#'   \item \strong{Annual WARM simulation:} Aggregate historical daily data to annual
#'     means (by water year if applicable), run wavelet analysis on the annual
#'     series, simulate candidate annual traces via \code{\link{wavelet_arima}},
#'     and subset traces using \code{\link{filter_warm_simulations}}.
#'   \item \strong{Daily disaggregation:} For each realization, call
#'     \code{\link{resample_weather_dates}} to generate a simulated sequence of
#'     historical analogue dates (in the internal 365-day calendar).
#'   \item \strong{Output construction:} Map resampled internal dates back to
#'     historical observation dates (\code{dateo}) and return simulated dates.
#' }
#'
#' @param obs_data Named list of data.frames (one per grid cell) with observed
#'   daily weather. Each data.frame must have one row per day corresponding to
#'   \code{obs_dates}, and columns for all \code{vars}.
#' @param obs_grid Data.frame describing grid cells. Must contain columns
#'   \code{xind}, \code{yind}, \code{x}, \code{y}. If \code{id} is missing, it is
#'   created as a sequential index.
#' @param obs_dates Vector of class \code{Date} corresponding to rows of each
#'   \code{obs_data[[i]]}. May include Feb 29; Feb 29 is removed internally along
#'   with the corresponding rows in \code{obs_data}.
#' @param vars Character vector of daily variables to simulate
#'   (e.g., \code{c("precip","temp")}).
#' @param n_years Integer number of years to simulate. If \code{NULL}, defaults
#'   to the number of complete historical simulation-years available after
#'   enforcing the 365-day calendar and restricting to the longest contiguous
#'   block of complete years.
#' @param start_year Integer first simulation year (calendar year if
#'   \code{year_start_month == 1}, otherwise first water year label).
#' @param year_start_month Integer in \code{1:12}. First month of the simulation year.
#'   Use \code{1} for calendar years; use \code{10} for Oct--Sep water years, etc.
#' @param n_realizations Integer number of realizations to generate.
#' @param warm_var Character name of the annual driver variable used in WARM.
#'   Must be present in \code{vars}.
#' @param warm_signif Numeric in (0,1). Wavelet significance level used to retain
#'   low-frequency components in WARM.
#' @param warm_pool_size Integer number of candidate annual traces to generate
#'   before filtering down to \code{n_realizations}.
#' @param annual_knn_n Integer number of historical analogue years considered in
#'   the annual KNN step of the disaggregation.
#' @param wet_q Numeric in (0,1). Wet-day threshold quantile used in the daily
#'   Markov chain classification.
#' @param extreme_q Numeric in (0,1). Extreme-day threshold quantile used in the
#'   daily Markov chain classification. Must be greater than \code{wet_q}.
#' @param dry_spell_factor Numeric length-12 vector of monthly dry-spell adjustment
#'   factors applied in the Markov-chain persistence logic.
#' @param wet_spell_factor Numeric length-12 vector of monthly wet-spell adjustment
#'   factors applied in the Markov-chain persistence logic.
#' @param out_dir Character directory for diagnostics and CSV outputs. The
#'   directory is created if it does not exist.
#' @param seed Optional integer seed for reproducibility. Used to initialize
#'   RNG streams for annual (WARM) and daily resampling components.
#' @param parallel Logical; if \code{TRUE}, run the daily disaggregation across
#'   realizations in parallel using a PSOCK cluster.
#' @param n_cores Optional integer number of cores for parallel execution. If
#'   \code{NULL} and \code{parallel = TRUE}, defaults to
#'   \code{max(1, parallel::detectCores() - 1)}.
#' @param verbose Logical; if \code{TRUE}, emits progress logs via \code{.log_info()}.
#'
#' @return A list with:
#' \describe{
#'   \item{resampled}{Tibble/data.frame with one column per realization. Each
#'     column contains historical observation dates (\code{dateo}) selected as
#'     analogues for each simulated day. Column names are \code{rlz_1}, \code{rlz_2}, ...}
#'   \item{dates}{Vector of class \code{Date} giving the simulated daily time axis
#'     (Feb 29 excluded).}
#' }
#'
#' @examples
#' \dontrun{
#' # obs_data: list of grid-cell data.frames with columns precip/temp
#' # obs_grid: data.frame with xind/yind/x/y
#' # obs_dates: Date vector aligned with rows of each obs_data[[i]]
#' out <- generate_weather(
#'   obs_data = obs_data,
#'   obs_grid = obs_grid,
#'   obs_dates = obs_dates,
#'   vars = c("precip", "temp"),
#'   start_year = 2020,
#'   n_years = 30,
#'   year_start_month = 10,
#'   n_realizations = 20,
#'   warm_var = "precip",
#'   warm_signif = 0.90,
#'   warm_pool_size = 5000,
#'   annual_knn_n = 120,
#'   wet_q = 0.30,
#'   extreme_q = 0.80,
#'   out_dir = tempdir(),
#'   seed = 123,
#'   parallel = TRUE,
#'   n_cores = 4,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @import dplyr
#' @import logger
#' @importFrom dplyr mutate
generate_weather <- function(
    obs_data = NULL,
    obs_grid = NULL,
    obs_dates = NULL,
    vars = NULL,
    n_years = NULL,
    start_year = 2020,
    year_start_month = 1,
    n_realizations = 5,
    warm_var = "precip",
    warm_signif = 0.90,
    warm_pool_size = 5000,
    annual_knn_n = 120,
    wet_q = 0.3,
    extreme_q = 0.8,
    dry_spell_factor = rep(1, 12),
    wet_spell_factor = rep(1, 12),
    out_dir = tempdir(),
    seed = NULL,
    parallel = FALSE,
    n_cores = NULL,
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
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  stopifnot(
    "obs_data must be a list" = is.list(obs_data),
    "obs_grid must be a data frame" = is.data.frame(obs_grid),
    "obs_dates must be a Date vector" = inherits(obs_dates, "Date"),
    "obs_data must be non-empty" = length(obs_data) >= 1L,
    "obs_dates length must match data rows" = length(obs_dates) == nrow(obs_data[[1]]),
    "year_start_month must be 1-12" = year_start_month %in% 1:12,
    "vars must be non-empty" = length(vars) > 0,
    "vars must exist in obs_data" = all(vars %in% names(obs_data[[1]])),
    "dry_spell_factor must have length 12" = length(dry_spell_factor) == 12,
    "wet_spell_factor must have length 12" = length(wet_spell_factor) == 12,
    "wet_q must be in (0,1)" = is.numeric(wet_q) && wet_q > 0 && wet_q < 1,
    "extreme_q must be in (0,1)" = is.numeric(extreme_q) && extreme_q > 0 && extreme_q < 1,
    "extreme_q > wet_q" = extreme_q > wet_q,
    "warm_var must be in vars" = warm_var %in% vars,
    "obs_grid must have required columns" = all(c("xind", "yind", "x", "y") %in% names(obs_grid))
  )

  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE/FALSE.", call. = FALSE)
  }

  if (length(unique(vapply(obs_data, nrow, integer(1)))) != 1L) {
    stop("All grid cells must have the same number of observations", call. = FALSE)
  }

  if (!("id" %in% names(obs_grid))) {
    obs_grid$id <- seq_len(nrow(obs_grid))
  }

  # ---------------------------------------------------------------------------
  # RNG management
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  warm_seed  <- sample.int(.Machine$integer.max, 1L)
  daily_seed <- sample.int(.Machine$integer.max, 1L)

  # ---------------------------------------------------------------------------
  # Parallel computing setup
  # ---------------------------------------------------------------------------
  if (isTRUE(parallel)) {
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = daily_seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    .log("Starting in parallel mode", tag = "INIT")
    .log("Number of cores: {n_cores}", tag = "INIT")
  } else {
    .log("Starting in sequential mode", tag = "INIT")
  }

  # ---------------------------------------------------------------------------
  # Initialization logging
  # ---------------------------------------------------------------------------
  n_grids <- length(obs_data)
  .log("Randomization seed: {seed}", tag = "INIT")
  .log("Climate variables: {paste(vars, collapse = ', ')}", tag = "INIT")
  .log("Total number of grids: {n_grids}", tag = "INIT")
  .log("Historical period: {obs_dates[1]} to {obs_dates[length(obs_dates)]}", tag = "INIT")

  # ---------------------------------------------------------------------------
  # Calendar normalization (enforce 365-day calendar)
  # ---------------------------------------------------------------------------
  leap_idx <- find_leap_days(obs_dates)
  if (!is.null(leap_idx) && length(leap_idx) > 0L) {
    obs_dates_old <- obs_dates
    obs_dates <- obs_dates[-leap_idx]

    obs_data <- lapply(obs_data, function(df) {
      if (!is.data.frame(df)) stop("Each element of obs_data must be a data.frame.", call. = FALSE)
      if (nrow(df) != length(obs_dates_old)) {
        stop("Each obs_data[[i]] must have nrow equal to length(obs_dates).", call. = FALSE)
      }
      df[-leap_idx, , drop = FALSE]
    })

    .log("Dropped {length(leap_idx)} row(s): 365-day calendar enforced", tag = "INIT")
  }

  # ---------------------------------------------------------------------------
  # Historical date table (internal time axis)
  # ---------------------------------------------------------------------------
  his_date  <- obs_dates
  his_wyear <- get_water_year(his_date, year_start_month)

  dates_d <- tibble::tibble(
    dateo = his_date,
    year  = as.integer(format(his_date, "%Y")),
    wyear = his_wyear,
    month = as.integer(format(his_date, "%m")),
    day   = as.integer(format(his_date, "%d"))
  )

  dates_d <- dates_d |>
    dplyr::mutate(
      date = if (year_start_month == 1) {
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

  wyear_idx <- match(dates_d$dateo, obs_dates)
  if (anyNA(wyear_idx)) stop("Internal error: dates_d$dateo did not match obs_dates.", call. = FALSE)

  .log("Using {length(full_wyears)} complete year(s) ({min(full_wyears)}-{max(full_wyears)})", tag = "INIT")

  # ---------------------------------------------------------------------------
  # Area-averaged daily and annual climate series
  # ---------------------------------------------------------------------------
  if (n_grids == 1L) {
    climate_d_aavg <- obs_data[[1]][wyear_idx, vars, drop = FALSE]
  } else {
    n_days <- length(wyear_idx)
    n_vars <- length(vars)

    climate_d_aavg_mat <- matrix(0, nrow = n_days, ncol = n_vars)
    colnames(climate_d_aavg_mat) <- vars

    for (i in seq_len(n_grids)) {
      grid_data <- as.matrix(obs_data[[i]][wyear_idx, vars, drop = FALSE])
      climate_d_aavg_mat <- climate_d_aavg_mat + grid_data
    }

    climate_d_aavg <- as.data.frame(climate_d_aavg_mat / n_grids)
  }

  climate_d_aavg$wyear <- dates_d$wyear

  climate_a_aavg <- climate_d_aavg |>
    dplyr::group_by(wyear) |>
    dplyr::summarize(dplyr::across(dplyr::all_of(vars), mean), .groups = "drop")

  # ---------------------------------------------------------------------------
  # Simulation date table (internal 365-day calendar)
  # ---------------------------------------------------------------------------
  if (is.null(n_years)) n_years <- length(unique(dates_d$wyear))
  n_years <- as.integer(n_years)

  sim_year_end   <- start_year + n_years
  sim_date_start <- as.Date(sprintf("%04d-%02d-01", start_year, year_start_month))
  sim_date_end   <- as.Date(sprintf("%04d-%02d-01", sim_year_end, year_start_month)) - 1

  sim_date_ini <- seq.Date(sim_date_start, sim_date_end, by = "day")
  sim_date_ini <- sim_date_ini[format(sim_date_ini, "%m-%d") != "02-29"]

  sim_dates_d <- tibble::tibble(
    dateo = sim_date_ini,
    year  = as.integer(format(sim_date_ini, "%Y")),
    wyear = get_water_year(sim_date_ini, year_start_month),
    month = as.integer(format(sim_date_ini, "%m")),
    day   = as.integer(format(sim_date_ini, "%d"))
  ) |>
    dplyr::mutate(
      date = if (year_start_month == 1) {
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
  warm_series_obs <- climate_a_aavg[[warm_var]]

  warm_power <- wavelet_spectral_analysis(
    variable = warm_series_obs,
    signif.level = warm_signif,
    period.lower.limit = 2,
    detrend = detrend,
    mode = "complete"
  )

  p <- plot_wavelet_spectra(
    variable = warm_series_obs,
    variable.year = climate_a_aavg$wyear,
    period = warm_power$gws_period,
    POWER = warm_power$power,
    GWS = warm_power$gws_unmasked,
    GWS_signif = warm_power$gws_signif_unmasked,
    coi = warm_power$coi,
    sigm = warm_power$sigm
  )

  tryCatch(
    ggplot2::ggsave(file.path(out_dir, "global_wavelet_power_spectrum.png"), p, width = 8, height = 5),
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
    sim.year.num = n_years,
    sim.num = warm_pool_size,
    seed = warm_seed,
    match.variance = TRUE
  )

  sim_annual_sub <- filter_warm_simulations(
    series.obs = warm_series_obs,
    series.sim = sim_annual,
    sample.num = n_realizations,
    seed = warm_seed + 1L,
    padding = TRUE,
    bounds = list(),
    relax.order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
    make.plots = TRUE,
    wavelet.pars = list(
      signif.level = warm_signif,
      noise.type = "red",
      period.lower.limit = 2,
      detrend = detrend
    ),
    verbose = verbose
  )

  tryCatch(
    {
      ggplot2::ggsave(file.path(out_dir, "warm_annual_series.png"), sim_annual_sub$plots[[1]], width = 8, height = 5)
      ggplot2::ggsave(file.path(out_dir, "warm_annual_statistics.png"), sim_annual_sub$plots[[2]], width = 8, height = 5)
      ggplot2::ggsave(file.path(out_dir, "warm_annual_wavelet.png"), sim_annual_sub$plots[[3]], width = 8, height = 5)
    },
    error = function(e) logger::log_warn("Failed to save warm plots: {e$message}")
  )

  # ---------------------------------------------------------------------------
  # Daily disaggregation (KNN + Markov chain) -> resampled historical dates
  # ---------------------------------------------------------------------------
  .log("Starting daily weather simulation using KNN + Markov Chain scheme", tag = "KNN")

  if (isTRUE(parallel)) {
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }

  resampled_ini <- foreach::foreach(n = seq_len(n_realizations)) %d% {
    resample_weather_dates(
      annual_prcp_sim = sim_annual_sub$sampled[, n],
      annual_prcp_obs = warm_series_obs,
      daily_prcp_obs  = climate_d_aavg$precip,
      daily_temp_obs  = climate_d_aavg$temp,
      sim_start_year  = start_year,
      realization_id  = n,
      n_sim_years     = n_years,
      obs_dates_table = dates_d,
      sim_dates_table = sim_dates_d,
      annual_knn_n    = annual_knn_n,
      year_start_month = year_start_month,
      wet_q           = wet_q,
      extreme_q       = extreme_q,
      dry_spell_factor = dry_spell_factor,
      wet_spell_factor = wet_spell_factor,
      seed = daily_seed + n
    )
  }

  # ---------------------------------------------------------------------------
  # Map internal simulated dates -> original historical observation dates (dateo)
  # ---------------------------------------------------------------------------
  resampled_dates <- vector("list", n_realizations)
  names(resampled_dates) <- paste0("rlz_", seq_len(n_realizations))

  for (i in seq_len(n_realizations)) {
    resampled_dates[[i]] <- dates_d$dateo[match(resampled_ini[[i]], dates_d$date)]
  }
  resampled_dates <- tibble::as_tibble(resampled_dates)

  # ---------------------------------------------------------------------------
  # Outputs
  # ---------------------------------------------------------------------------
  tryCatch(
    utils::write.csv(sim_dates_d$date, file.path(out_dir, "sim_dates.csv"), row.names = FALSE),
    error = function(e) logger::log_error("Failed to write sim_dates.csv: {e$message}")
  )

  tryCatch(
    utils::write.csv(resampled_dates, file.path(out_dir, "resampled_dates.csv"), row.names = FALSE),
    error = function(e) logger::log_error("Failed to write resampled_dates.csv: {e$message}")
  )

  .log("Results saved to: `{out_dir}`", tag = "DONE")
  .log("Simulation complete. Elapsed time: {Sys.time() - start_time} secs", tag = "DONE")

  list(resampled = resampled_dates, dates = sim_dates_d$dateo)
}

