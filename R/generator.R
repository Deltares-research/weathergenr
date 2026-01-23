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
#'     series, simulate candidate annual traces via \code{\link{simulate_warm}},
#'     and subset traces using \code{\link{filter_warm_pool}}.
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
#' @param warm_filter_bounds Named list of filtering thresholds and relaxation controls
#'   forwarded to \code{\link{filter_warm_pool}} as \code{filter_bounds}. Any entry
#'   overrides internal defaults. Uses snake_case keys (e.g. \code{tail_low_p}).
#' @param relax_priority Character vector giving the relaxation priority
#'   for WARM filtering, forwarded to \code{\link{filter_warm_pool}} as \code{relax_order}.
#'   Must contain each of \code{c("mean","sd","tail_low","tail_high","wavelet")} exactly once.
#'   Filters are relaxed iteratively by loosening the currently most restrictive criterion
#'   (lowest pass rate), subject to this priority ordering.
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
#' @param verbose Logical; if \code{TRUE}, emits progress logs via \code{.log()}.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
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
#' @importFrom dplyr mutate
#' @importFrom foreach %dopar%
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
  warm_filter_bounds = list(),
  relax_priority = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
  annual_knn_n = 120,
  wet_q = 0.3,
  extreme_q = 0.8,
  dry_spell_factor = rep(1, 12),
  wet_spell_factor = rep(1, 12),
  out_dir = tempdir(),
  seed = NULL,
  parallel = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  save_plots = TRUE
) {


  start_time <- Sys.time()
  verbose <- isTRUE(verbose)
  save_plots <- isTRUE(save_plots)

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
    "All grid cells must have the same number of rows" = length(unique(vapply(obs_data, nrow, integer(1)))) == 1L,
    "year_start_month must be 1-12" = year_start_month %in% 1:12,
    "vars must be non-empty" = length(vars) > 0,
    "vars must exist in obs_data" = all(vars %in% names(obs_data[[1]])),
    "dry_spell_factor must have length 12" = length(dry_spell_factor) == 12,
    "wet_spell_factor must have length 12" = length(wet_spell_factor) == 12,
    "wet_q must be in (0,1)" = is.numeric(wet_q) && wet_q > 0 && wet_q < 1,
    "extreme_q must be in (0,1)" = is.numeric(extreme_q) && extreme_q > 0 && extreme_q < 1,
    "extreme_q must be greater than wet_q" = extreme_q > wet_q,
    "warm_var must be in vars" = warm_var %in% vars,
    "obs_grid must have required columns" = all(c("xind", "yind", "x", "y") %in% names(obs_grid)),
    "verbose must be TRUE or FALSE" = is.logical(verbose) && length(verbose) == 1L,
    "warm_filter_bounds must be a list" = is.list(warm_filter_bounds),
    "relax_priority must be a character vector" = is.character(relax_priority)
  )


  if (!("id" %in% names(obs_grid))) obs_grid$id <- seq_len(nrow(obs_grid))

  # ---------------------------------------------------------------------------
  # CONSTANTS -- no user-level change needed for now
  # ---------------------------------------------------------------------------

  ## GGPLOT constants
  pl_width <- 8
  pl_height <- 5

  # WARM constants
  DETREND <- TRUE
  MIN_PERIOD <- 2
  WARM_FILTER <- "la8"


  # ---------------------------------------------------------------------------
  # RNG management
  # ---------------------------------------------------------------------------

  # save and restore state
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_rng_state <- .Random.seed
      has_rng_state <- TRUE
    } else {
      has_rng_state <- FALSE
    }
    on.exit({
      if (has_rng_state) .Random.seed <<- old_rng_state
    }, add = TRUE)
    set.seed(seed)
  }

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

    .log("Starting in parallel mode", tag = "INIT", verbose = verbose)
    .log("Number of cores: {n_cores}", tag = "INIT", verbose = verbose)
  } else {
    .log("Starting in sequential mode", tag = "INIT", verbose = verbose)
  }

  # ---------------------------------------------------------------------------
  # Initial Logging
  # ---------------------------------------------------------------------------

  n_grids <- length(obs_data)

  # General info
  .log("Historical period: {obs_dates[1]} to {obs_dates[length(obs_dates)]}", tag = "INIT", verbose = verbose)
  .log(paste0("Variables: ", paste(as.character(vars), collapse = ", ")), tag = "INIT", verbose = verbose)
  .log("Grid cells: {n_grids}", tag = "INIT", verbose = verbose)

  # Calendar / water-year logic
  .log("Year logic: {if (as.integer(year_start_month) == 1L) 'calendar year (Jan-Dec)' else paste0('water year (start month = ', as.integer(year_start_month), ')')}",
       tag = "INIT", verbose = verbose)

  # Output directory
  .log("Output directory: {normalizePath(out_dir, winslash = '/', mustWork = FALSE)}",
       tag = "INIT", verbose = verbose)

  # Seed information
  if (!is.null(seed)) {
    .log("Seed: {seed}", tag = "INIT", verbose = verbose)
  } else {
    .log("Seed: not set (non-reproducible)", tag = "INIT", verbose = verbose)
  }

  # ---------------------------------------------------------------------------
  # Preprocessing
  # ---------------------------------------------------------------------------

  # 4a. Normalize calendar (remove leap days)
  cal_norm <- normalize_calendar(obs_dates, obs_data, verbose = verbose)
  obs_dates <- cal_norm$dates
  obs_data  <- cal_norm$data

  # 4b. Build historical date table
  hist_dates <- build_historical_dates(obs_dates, year_start_month, verbose = verbose)
  dates_d      <- hist_dates$dates_df
  wyear_idx    <- hist_dates$wyear_idx

  # 4c. Compute area averages
  area_avg <- compute_area_averages(obs_data, wyear_idx, dates_d$wyear, vars)
  climate_d_aavg <- area_avg$daily
  climate_a_aavg <- area_avg$annual

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
    wyear = compute_water_year(sim_date_ini, year_start_month),
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
  # Annual time-series generation (WARM) - Hybrid CWT/MODWT approach
  # ---------------------------------------------------------------------------
  .log("Running annual WARM simulation", tag = "WARM", verbose = verbose)

  # Mean-center (anomaly by mean only)
  warm_series_obs <- climate_a_aavg[[warm_var]]
  warm_series_obs_mean <- mean(warm_series_obs)
  warm_series_obs_anom <- warm_series_obs - warm_series_obs_mean

  # Hybrid analysis: CWT for visualization, MODWT MRA for additive components
  warm_analysis <- analyze_wavelet_additive(
    series = warm_series_obs_anom,
    signif = warm_signif,
    noise = "red",
    min_period = MIN_PERIOD,
    detrend = DETREND,
    filter = WARM_FILTER,
    n_levels = NULL,
    include_smooth = TRUE,
    cwt_mode = "complete",
    diagnostics = FALSE
  )

  if(save_plots) {

    # Use CWT results for power spectrum visualization
    p <- plot_wavelet_power(
      series = warm_series_obs,
      time = climate_a_aavg$wyear,
      period = warm_analysis$cwt$period,
      power = warm_analysis$cwt$power,
      gws = warm_analysis$cwt$gws_unmasked,
      gws_signif = warm_analysis$cwt$gws_signif_unmasked,
      coi = warm_analysis$cwt$coi,
      signif_mask = warm_analysis$cwt$sigm
    )

    tryCatch(
      ggsave(file.path(out_dir, "obs_power_spectra.png"), p, width = 8, height = 5),
      error = function(e) .log("Failed to save wavelet spectra plot: {e$message}", level = "warn", verbose = verbose)
    )

  }

  # Log significant periodicities from CWT analysis
  if (warm_analysis$has_significance) {
    .log("Significant periodicities (CWT): {length(warm_analysis$cwt_signif_periods)}", tag = "WARM", verbose = verbose)
    .log("Periods (years): {paste(round(warm_analysis$cwt_signif_periods, 2), collapse = ', ')}", tag = "WARM", verbose = verbose)
  } else {
    .log("No significant low-frequency periodicities detected", tag = "WARM", verbose = verbose)
  }

  # Log MODWT MRA component information
  .log("Wavelet (MODWT): {ncol(warm_analysis$components)} components ({paste(warm_analysis$component_names, collapse = ', ')})",
    tag = "WARM", verbose = verbose)

  # Use MODWT MRA additive components for WARM simulation
  sim_annual_anom <- simulate_warm(
    components = warm_analysis$components,
    n = n_years,
    n_sim = warm_pool_size,
    seed = warm_seed,
    match_variance = TRUE,
    verbose = verbose)

  # Add the mean back to get absolute annual precipitation
  sim_annual <- sim_annual_anom + warm_series_obs_mean

  sim_annual_sub <- filter_warm_pool(
    obs_series = warm_series_obs,
    sim_series = sim_annual,
    n_select = n_realizations,
    seed = warm_seed + 1L,
    filter_bounds = warm_filter_bounds,
    relax_order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
    make_plots = TRUE,
    wavelet_args = list(
      signif_level = warm_signif,
      noise_type = "red",
      period_lower_limit = 2,
      detrend = DETREND),
    verbose = verbose,
    parallel = parallel,
    n_cores  = n_cores)

  if(save_plots) {
    tryCatch(
      {
        ggsave(file.path(out_dir, "warm_annual_precip.png"),
               sim_annual_sub$plots[[1]], width = pl_width, height = pl_height)
        ggsave(file.path(out_dir, "warm_annual_stats.png"),
               sim_annual_sub$plots[[2]], width = pl_width, height = pl_height)
        ggsave(file.path(out_dir, "warm_annual_wavelet.png"),
               sim_annual_sub$plots[[3]], width = pl_width, height = pl_height)
      },
      error = function(e) .log("Failed to save warm filtering plots: {e$message}", level = "warn")
    )
  }

  # ---------------------------------------------------------------------------
  # Daily disaggregation (KNN + Markov chain) -> resampled historical dates
  # ---------------------------------------------------------------------------
  .log("Running daily KNN + Markov Chain resampling", tag = "RESAMPLE", verbose = verbose)

  if (isTRUE(parallel)) {

    # Parallel: foreach %dopar%
    .log("Processing {n_realizations} realizations...", tag = "RESAMPLE", verbose = verbose)

    resampled_ini <- foreach::foreach(n = seq_len(n_realizations)) %dopar% {
      resample_weather_dates(
        sim_annual_precip = sim_annual_sub$selected[, n],
        obs_annual_precip = warm_series_obs,
        obs_daily_precip  = climate_d_aavg$precip,
        obs_daily_temp  = climate_d_aavg$temp,
        year_start  = start_year,
        n_years     = n_years,
        obs_dates_df = dates_d,
        sim_dates_df = sim_dates_d,
        annual_knn_n    = annual_knn_n,
        year_start_month = year_start_month,
        wet_q           = wet_q,
        extreme_q       = extreme_q,
        dry_spell_factor = dry_spell_factor,
        wet_spell_factor = wet_spell_factor,
        seed = daily_seed + n)
    }

  } else {

    # Sequential: plain for loop with progress
    resampled_ini <- vector("list", n_realizations)
    progress_parts <- character(n_realizations)

    for (n in seq_len(n_realizations)) {

      progress_parts[n] <- paste0(n, "/", n_realizations)

      resampled_ini[[n]] <-  resample_weather_dates(
        sim_annual_precip = sim_annual_sub$selected[, n],
        obs_annual_precip = warm_series_obs,
        obs_daily_precip  = climate_d_aavg$precip,
        obs_daily_temp  = climate_d_aavg$temp,
        year_start  = start_year,
        n_years     = n_years,
        obs_dates_df = dates_d,
        sim_dates_df = sim_dates_d,
        annual_knn_n    = annual_knn_n,
        year_start_month = year_start_month,
        wet_q  = wet_q,
        extreme_q   = extreme_q,
        dry_spell_factor = dry_spell_factor,
        wet_spell_factor = wet_spell_factor,
        seed = daily_seed + n)
    }

    progress_str <- paste(progress_parts, collapse = "..")
    .log("Processing realization: {progress_str}", tag = "RESAMPLE", verbose = verbose)
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
    utils::write.csv(data.frame(date = sim_dates_d$date), file.path(out_dir, "sim_dates.csv"), row.names = FALSE),
    error = function(e) {
      msg <- paste0("Failed to write sim_dates.csv: ", conditionMessage(e))
      .log(msg, level = "error", verbose = verbose)
    }
  )

  tryCatch(
    utils::write.csv(resampled_dates, file.path(out_dir, "resampled_dates.csv"), row.names = FALSE),
    error = function(e) {
      msg <- paste0("Failed to write resampled_dates.csv: ", conditionMessage(e))
      .log(msg, level = "error", verbose = verbose)
    }
  )

  .log("Output: `{out_dir}`", tag = "COMPLETE", verbose = verbose)
  .log("Elapsed time: {format_elapsed(start_time)}", tag = "COMPLETE", verbose = verbose)

  list(resampled = resampled_dates, dates = sim_dates_d$dateo)
}



#' Run weathergenr end-to-end (generate and evaluate)
#'
#' @description
#' Top-level convenience wrapper that runs:
#' \code{generate_weather()}
#' \code{prepare_evaluation_data()}
#' \code{evaluate_weather_generator()}
#'
#' All execution and evaluation settings are taken from \code{config}.
#' Optional logging writes console output to a timestamped file in \code{out_dir}
#' while continuing to display output on the console.
#'
#' @param obs_data Observed data (e.g. \code{ncdata$data}).
#' @param obs_grid Observed grid metadata (e.g. \code{ncdata$grid}).
#' @param obs_dates Observed dates (e.g. \code{ncdata$date}).
#' @param out_dir Character. Output directory.
#' @param config List. Full simulation/evaluation configuration.
#' @param eval_max_grids Integer. Maximum number of grids to evaluate.
#' @param log_messages Logical. If TRUE, save console output to
#'   \code{log_YYYYMMDD_HHMMSS.txt} in \code{out_dir}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{gen_output}: output of \code{generate_weather()}
#'   \item \code{evaluation}: output of \code{evaluate_weather_generator()}
#'   \item \code{log_path}: path to the log file (or NULL if \code{log_messages=FALSE})
#' }
#'
#' @export
run_weather_generator <- function(
    obs_data,
    obs_grid,
    obs_dates,
    out_dir,
    config,
    eval_max_grids = 25L,
    log_messages = FALSE
) {
  if (!is.list(config)) {
    stop("'config' must be a list.", call. = FALSE)
  }
  if (is.null(config$vars) || length(config$vars) == 0) {
    stop("'config$vars' must be provided.", call. = FALSE)
  }
  if (is.null(config$n_realizations)) {
    stop("'config$n_realizations' must be provided.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Optional console logging (console + file)
  # ---------------------------------------------------------------------------

  .strip_ansi <- function(x) gsub("\033\\[[0-9;]*[A-Za-z]", "", x)

  log_path <- NULL
  con <- NULL
  msg_con <- NULL
  msg_file <- NULL

  if (isTRUE(log_messages)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (file.access(out_dir, 2) != 0) {
      stop(
        "Cannot write to out_dir: ",
        normalizePath(out_dir, winslash = "/", mustWork = FALSE),
        call. = FALSE
      )
    }

    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_path <- file.path(out_dir, paste0("log_", ts, ".txt"))
    msg_file <- file.path(out_dir, paste0("log_", ts, "_messages.txt"))

    # Main log connection (regular output via sink)
    con <- file(log_path, open = "wt", encoding = "UTF-8")
    writeLines(
      c(
        "weathergenr run log",
        paste0("started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0("out_dir: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
        ""
      ),
      con = con,
      useBytes = TRUE
    )
    flush(con)

    # Mirror regular output (cat/print) to file while keeping console output visible
    sink(con, split = TRUE)

    # Messages/warnings are captured separately via handlers (Windows limitation)
    msg_con <- file(msg_file, open = "wt", encoding = "UTF-8")
    writeLines(
      c(
        "MESSAGES/WARNINGS (captured separately)",
        ""
      ),
      con = msg_con,
      useBytes = TRUE
    )
    flush(msg_con)

    on.exit({
      # Unwind sink stack (regular output only)
      while (sink.number() > 0) sink(NULL)

      if (!is.null(con) && isOpen(con)) close(con)
      if (!is.null(msg_con) && isOpen(msg_con)) close(msg_con)

      # Sanitize and merge messages into main log, then remove message temp file
      if (!is.null(log_path) && file.exists(log_path)) {
        main <- readLines(log_path, warn = FALSE, encoding = "UTF-8")
        main <- .strip_ansi(main)
        main <- gsub("\r$", "", main)

        if (!is.null(msg_file) && file.exists(msg_file)) {
          msgs <- readLines(msg_file, warn = FALSE, encoding = "UTF-8")
          msgs <- .strip_ansi(msgs)
          msgs <- gsub("\r$", "", msgs)

          out <- c(main, "----- MESSAGES/WARNINGS -----", msgs)
          writeLines(out, log_path, useBytes = TRUE)

          suppressWarnings(file.remove(msg_file))
        } else {
          writeLines(main, log_path, useBytes = TRUE)
        }
      }
    }, add = TRUE)
  }

  # ---------------------------------------------------------------------------
  # Run body (wrapped so we can tee messages/warnings without modifying internals)
  # ---------------------------------------------------------------------------

  run_body <- function() {
    # -------------------------------------------------------------------------
    # Derive grid IDs
    # -------------------------------------------------------------------------
    derive_grid_ids <- function(grid) {
      if (is.list(grid) && !is.null(grid$id)) return(grid$id)
      if (is.data.frame(grid) && "id" %in% names(grid)) return(grid$id)

      n <- NA_integer_
      if (is.data.frame(grid)) n <- nrow(grid)
      if (is.matrix(grid)) n <- nrow(grid)
      if (is.list(grid) && !is.null(grid$lon) && !is.null(grid$lat)) {
        if (is.numeric(grid$lon) && is.numeric(grid$lat)) {
          n <- length(grid$lon)
        }
      }

      if (is.na(n) || n <= 0L) {
        stop(
          "Unable to derive grid IDs from 'obs_grid'. ",
          "Provide 'obs_grid$id' or a data.frame with column 'id'.",
          call. = FALSE
        )
      }
      seq_len(n)
    }

    grid_ids <- derive_grid_ids(obs_grid)

    # -------------------------------------------------------------------------
    # Step 1: Generate synthetic weather
    # -------------------------------------------------------------------------
    gen_output <- generate_weather(
      obs_data           = obs_data,
      obs_grid           = obs_grid,
      obs_dates          = obs_dates,
      vars               = config$vars,
      n_years            = config$n_years,
      start_year         = config$start_year,
      year_start_month   = config$year_start_month,
      n_realizations     = config$n_realizations,
      warm_var           = config$warm_var,
      warm_signif        = config$warm_signif,
      warm_pool_size     = config$warm_pool_size,
      warm_filter_bounds = config$warm_filter_bounds,
      annual_knn_n       = config$annual_knn_n,
      wet_q              = config$wet_q,
      extreme_q          = config$extreme_q,
      dry_spell_factor   = config$dry_spell_factor,
      wet_spell_factor   = config$wet_spell_factor,
      out_dir            = out_dir,
      parallel           = config$parallel,
      n_cores            = config$n_cores,
      seed               = config$seed,
      verbose            = config$verbose,
      save_plots         = config$save_plots
    )

    # -------------------------------------------------------------------------
    # Step 2: Prepare evaluation data (internal)
    # -------------------------------------------------------------------------
    eval_data <- prepare_evaluation_data(
      gen_output = gen_output,
      obs_data   = obs_data,
      obs_dates  = obs_dates,
      grid_ids   = grid_ids,
      variables  = config$vars,
      verbose    = isTRUE(config$verbose)
    )

    # -------------------------------------------------------------------------
    # Step 3: Evaluate generator performance
    # -------------------------------------------------------------------------
    evaluation <- evaluate_weather_generator(
      daily_sim       = eval_data$sim_data,
      daily_obs       = eval_data$obs_data,
      vars            = config$vars,
      variable_labels = NULL,
      n_realizations  = config$n_realizations,
      eval_max_grids  = eval_max_grids,
      wet_q           = config$wet_q,
      extreme_q       = config$extreme_q,
      output_dir      = out_dir,
      save_plots      = isTRUE(config$save_plots),
      seed            = config$seed
    )

    list(
      gen_output = gen_output,
      evaluation = evaluation
    )
  }

  if (!isTRUE(log_messages)) {
    res <- run_body()
    res$log_path <- NULL
    return(res)
  }

  # With logging: tee messages/warnings to msg_file while keeping console output.
  res <- withCallingHandlers(
    run_body(),

    message = function(m) {
      # Let R print the message normally (do NOT call message() here)
      txt <- conditionMessage(m)

      # Write a clean single line to msg_file
      if (!is.null(msg_con) && isOpen(msg_con)) {
        txt2 <- sub("[\r\n]+$", "", txt)          # remove trailing newline(s) only
        cat(.strip_ansi(txt2), "\n", file = msg_con)
        flush(msg_con)
      }

      # Do NOT muffle -> message continues to console once
      NULL
    },

    warning = function(w) {
      # Let R print the warning normally (do NOT call warning() here)
      txt <- conditionMessage(w)

      if (!is.null(msg_con) && isOpen(msg_con)) {
        txt2 <- sub("[\r\n]+$", "", txt)
        cat("WARNING: ", .strip_ansi(txt2), "\n", file = msg_con, sep = "")
        flush(msg_con)
      }

      # Do NOT muffle -> warning continues to console once
      NULL
    }
  )


  res$log_path <- log_path
  res
}
