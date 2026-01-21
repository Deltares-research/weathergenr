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
#' @param warm_filter_relax_order Character vector giving the relaxation priority
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
  warm_filter_relax_order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
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
    "warm_filter_relax_order must be a character vector" = is.character(warm_filter_relax_order)
  )


  if (!("id" %in% names(obs_grid))) obs_grid$id <- seq_len(nrow(obs_grid))

  # ---------------------------------------------------------------------------
  # CONSTANTS
  # ---------------------------------------------------------------------------

  ## GGPLOT constants
  pl_width <- 8
  pl_height <- 5

  # WARM constants
  DETREND <- TRUE
  MIN_PERIOD <- 2

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
  if (!is.null(seed)) {
    .log("Seed: {seed}", tag = "INIT", verbose = verbose)
  } else {
    .log("Seed: not set (non-reproducible)", tag = "INIT", verbose = verbose)
  }
  .log(paste0("Variables: ", paste(as.character(vars), collapse = ", ")), tag = "INIT", verbose = verbose)
  .log("Grid cells: {n_grids}", tag = "INIT", verbose = verbose)
  .log("Historical period: {obs_dates[1]} to {obs_dates[length(obs_dates)]}", tag = "INIT", verbose = verbose)

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
  # Annual time-series generation (WARM)
  # ---------------------------------------------------------------------------
  .log("Running annual WARM simulation", tag = "WARM", verbose = verbose)


  warm_series_obs <- climate_a_aavg[[warm_var]]

  warm_power <- analyze_wavelet_spectrum(
    series = warm_series_obs,
    signif = warm_signif,
    noise = "red",
    min_period = MIN_PERIOD,
    detrend = DETREND,
    mode = "complete")

  p <- plot_wavelet_power(
    series = warm_series_obs,
    time = climate_a_aavg$wyear,
    period = warm_power$period,
    power = warm_power$power,
    gws = warm_power$gws_unmasked,
    gws_signif = warm_power$gws_signif_unmasked,
    coi = warm_power$coi,
    signif_mask = warm_power$sigm)

  tryCatch(
    ggsave(file.path(out_dir, "global_wavelet_power_spectrum.png"), p, width = 8, height = 5),
    error = function(e) .log("Failed to save wavelet spectrum plot: {e$message}", level = "warn", verbose = verbose)
  )

  if (any(!is.na(warm_power$signif_periods))) {
    .log("Significant low-frequency components: {length(warm_power$comps_names) - 1}", tag = "WARM", verbose = verbose)
    .log("Periodicity (years): {paste(round(warm_power$period[warm_power$signif_periods], 2), collapse = ', ')}",
      tag = "WARM", verbose = verbose)
  } else {
    .log("No significant low-frequency periodicities detected", tag = "WARM", verbose = verbose)
  }

  sim_annual <- simulate_warm(
    components = warm_power$comps,
    n = n_years,
    n_sim = warm_pool_size,
    seed = warm_seed,
    match_variance = TRUE,
    verbose = verbose
  )

  sim_annual_sub <- filter_warm_pool(
    obs_series = warm_series_obs,
    sim_series = sim_annual,
    n_select = n_realizations,
    seed = warm_seed + 1L,
    pad_periods = TRUE,
    filter_bounds = warm_filter_bounds,
    relax_order = warm_filter_relax_order,
    make_plots = TRUE,
    wavelet_args = list(
      signif_level = warm_signif,
      noise_type = "red",
      period_lower_limit = 2,
      detrend = DETREND
    ),
    verbose = verbose
  )

  tryCatch(
    {
      ggsave(file.path(out_dir, "warm_annual_series.png"),
             sim_annual_sub$plots[[1]], width = pl_width, height = pl_height)
      ggsave(file.path(out_dir, "warm_annual_statistics.png"),
             sim_annual_sub$plots[[2]], width = pl_width, height = pl_height)
      ggsave(file.path(out_dir, "warm_annual_wavelet.png"),
             sim_annual_sub$plots[[3]], width = pl_width, height = pl_height)
    },
    error = function(e) .log("Failed to save warm plots: {e$message}", level = "warn")
  )

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
