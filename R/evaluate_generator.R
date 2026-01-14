#' Evaluate stochastic weather generator performance (multi-grid, multi-realization)
#'
#' @description
#' Comprehensive diagnostic evaluation comparing **synthetic daily simulations** against
#' **historical daily observations** across multiple grid cells. The workflow:
#' \itemize{
#'   \item standardizes obs/sim periods to **full 365-day years** (leap days removed) and
#'         matches lengths via random contiguous windows,
#'   \item computes summary statistics (mean/sd/skewness), wet/dry day counts, wet/dry spell
#'         lengths, and correlation diagnostics,
#'   \item generates diagnostic plots and an overall fit summary across realizations.
#' }
#'
#' Variable naming is aligned with the package API:
#' precipitation is expected as \code{precip}; temperature variables can be evaluated as
#' \code{temp}, \code{temp_min}, \code{temp_max}, etc.
#'
#' @param sim_data List of synthetic realizations. Each element is a list of data.frames
#'   (one per grid cell). Each grid data.frame must include a \code{date} column and
#'   the variables listed in \code{vars}.
#' @param obs_data List of observed grid data.frames (one per grid cell). Each must include
#'   a \code{date} column and the variables listed in \code{vars}.
#' @param vars Character vector of variables to evaluate (e.g., \code{c("precip", "temp")}).
#'   Must include \code{"precip"} to enable wet/dry day and spell diagnostics.
#' @param n_realizations Integer. Number of realizations in \code{sim_data}. Must equal \code{length(sim_data)}.
#' @param wet_q Numeric in (0, 1). Monthly quantile threshold to define wet days from positive
#'   precipitation values (default = 0.2).
#' @param extreme_q Numeric in (0, 1). Monthly quantile threshold to define extremely wet days
#'   from positive precipitation values (default = 0.8). Must be > \code{wet_q}.
#' @param out_dir Character. Directory to save plots. If \code{NULL}, plots are not saved.
#' @param save_plot Logical. If \code{TRUE}, save plots to \code{out_dir}.
#' @param show_title Logical. If \code{TRUE}, include plot titles.
#' @param max_grid Integer. Maximum number of grid cells to evaluate (default = 25).
#'   If \code{length(obs_data) > max_grid}, a random subset is used to control memory.
#' @param seed Optional integer. If provided, ensures reproducible grid subsampling and
#'   year-window selection.
#'
#' @return A named list of \code{ggplot2} objects with class \code{"weather_assessment"}.
#'   The object includes attributes:
#'   \itemize{
#'     \item \code{fit_summary}: per-realization metrics and overall ranking
#'     \item \code{metadata}: run configuration summary
#'   }
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom e1071 skewness
#' @importFrom stats cor sd quantile
#' @importFrom utils combn
#' @export
evaluate_weather_generator <- function(
    sim_data = NULL,
    obs_data = NULL,
    vars = NULL,
    n_realizations = NULL,
    wet_q = 0.2,
    extreme_q = 0.8,
    out_dir = NULL,
    save_plot = TRUE,
    show_title = TRUE,
    max_grid = 25,
    seed = NULL
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  validate_weather_generator_inputs(
    sim_data = sim_data,
    obs_data = obs_data,
    vars = vars,
    n_realizations = n_realizations,
    wet_q = wet_q,
    extreme_q = extreme_q,
    max_grid = max_grid,
    seed = seed
  )

  # ---------------------------------------------------------------------------
  # RNG state management (only if seed provided)
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }
    on.exit({ if (has_seed) .Random.seed <<- old_seed }, add = TRUE)
    set.seed(as.integer(seed))
  }

  # ---------------------------------------------------------------------------
  # Logging helper
  # ---------------------------------------------------------------------------
  .log <- function(level = c("info", "warn"), msg) {
    level <- match.arg(level)
    if (!requireNamespace("logger", quietly = TRUE)) return(invisible(NULL))
    if (level == "info") logger::log_info(msg) else logger::log_warn(msg)
    invisible(NULL)
  }

  # ---------------------------------------------------------------------------
  # Setup
  # ---------------------------------------------------------------------------
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    save_plot <- FALSE
  }

  options(dplyr.summarise.inform = FALSE)

  plot_config <- list(
    subtitle = "Observed vs simulated distributions (all realizations)",
    alpha = 0.4,
    colors = c(Observed = "blue", Simulated = "gray40"),
    theme = theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 10)
      )
  )

  .log("info", sprintf("[ASSESS] Start | grids=%d | realizations=%d | vars=%s",
                       length(obs_data), n_realizations, paste(vars, collapse = ", ")))
  .log("info", sprintf("[ASSESS] Parameters | wet_q=%.3f | extreme_q=%.3f | max_grid=%d",
                       wet_q, extreme_q, max_grid))

  # ---------------------------------------------------------------------------
  # Grid subsampling (memory control)
  # ---------------------------------------------------------------------------
  n_grid_org <- length(obs_data)
  if (n_grid_org > max_grid) {
    sel <- sort(sample.int(n_grid_org, size = max_grid))
    obs_data <- obs_data[sel]
    sim_data <- lapply(sim_data, function(rlz) rlz[sel])
    .log("warn", sprintf("[ASSESS] Grid count reduced from %d to %d for memory control",
                         n_grid_org, length(obs_data)))
  }
  n_grid <- length(obs_data)

  # ---------------------------------------------------------------------------
  # Standardize periods: remove leap days + match full-year windows
  # ---------------------------------------------------------------------------
  .log("info", "[ASSESS] Standardizing obs/sim periods to full years (365-day) and equal length")
  std <- standardize_obs_sim_periods(
    obs_data = obs_data,
    sim_data = sim_data,
    n_realizations = n_realizations
  )
  obs_data <- std$obs_data
  sim_data <- std$sim_data

  .log("info", sprintf("[ASSESS] Standardized period | Obs=%d-%d | Sim=%d-%d | years=%d",
                       std$obs_year_start, std$obs_year_end,
                       std$sim_year_start, std$sim_year_end,
                       std$n_year))

  # ---------------------------------------------------------------------------
  # Process observed
  # ---------------------------------------------------------------------------
  .log("info", "[ASSESS] Processing observed data")
  obs_res <- process_observed_data(
    obs_data = obs_data,
    vars = vars,
    n_grid = n_grid,
    wet_q = wet_q,
    extreme_q = extreme_q
  )

  # ---------------------------------------------------------------------------
  # Process simulated
  # ---------------------------------------------------------------------------
  .log("info", sprintf("[ASSESS] Processing simulated data (%d realizations)", n_realizations))
  sim_res <- process_simulated_data(
    sim_data = sim_data,
    n_realizations = n_realizations,
    vars = vars,
    mc_thresholds = obs_res$mc_thresholds
  )

  # ---------------------------------------------------------------------------
  # Merge for plotting
  # ---------------------------------------------------------------------------
  .log("info", "[ASSESS] Preparing diagnostic data for plotting")
  plot_dat <- prepare_plot_data(
    obs_results = obs_res,
    sim_results = sim_res,
    vars = vars
  )

  # ---------------------------------------------------------------------------
  # Plots
  # NOTE: This expects your existing plotting factory to exist as-is.
  # ---------------------------------------------------------------------------
  .log("info", "[ASSESS] Generating diagnostic plots")
  plots <- create_all_diagnostic_plots(
    plot_data = plot_dat,
    plot_config = plot_config,
    vars = vars,
    show_title = show_title,
    save_plot = save_plot,
    out_dir = out_dir
  )

  .log("info", sprintf("[ASSESS] Generated %d diagnostic plots", length(plots)))
  if (isTRUE(save_plot)) .log("info", sprintf("[ASSESS] Plots saved to: %s", out_dir))

  # ---------------------------------------------------------------------------
  # Fit summary
  # ---------------------------------------------------------------------------
  .log("info", "[ASSESS] Computing fit metrics for all realizations")
  fit_summary <- compute_realization_fit_metrics(
    obs_results = obs_res,
    sim_results = sim_res,
    vars = vars
  )

  .log("info", "[ASSESS] Displaying fit assessment summary")
  display_fit_summary_table(fit_summary, vars)

  .log("info", "[ASSESS] Completed")

  structure(
    plots,
    class = c("weather_assessment", "list"),
    fit_summary = fit_summary,
    metadata = list(
      n_grid = n_grid,
      n_realizations = n_realizations,
      vars = vars,
      assessment_date = Sys.Date()
    )
  )
}

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

#' Validate inputs for evaluate_weather_generator()
#' @keywords internal
validate_weather_generator_inputs <- function(sim_data, obs_data, vars, n_realizations,
                                              wet_q, extreme_q, max_grid, seed) {

  if (is.null(sim_data)) stop("'sim_data' must not be NULL", call. = FALSE)
  if (is.null(obs_data)) stop("'obs_data' must not be NULL", call. = FALSE)
  if (is.null(vars)) stop("'vars' must not be NULL", call. = FALSE)
  if (is.null(n_realizations)) stop("'n_realizations' must not be NULL", call. = FALSE)

  if (!is.list(sim_data)) stop("'sim_data' must be a list", call. = FALSE)
  if (!is.list(obs_data)) stop("'obs_data' must be a list", call. = FALSE)
  if (!is.character(vars) || length(vars) < 1L) stop("'vars' must be a character vector", call. = FALSE)

  if (!is.numeric(n_realizations) || length(n_realizations) != 1L || !is.finite(n_realizations) ||
      n_realizations < 1 || (n_realizations %% 1) != 0) {
    stop("'n_realizations' must be a positive integer", call. = FALSE)
  }
  n_realizations <- as.integer(n_realizations)

  if (length(sim_data) != n_realizations) {
    stop("Length of 'sim_data' must equal 'n_realizations'", call. = FALSE)
  }
  if (length(obs_data) < 1L) stop("'obs_data' must contain at least one grid cell", call. = FALSE)

  # Required for wet/dry diagnostics
  if (!"precip" %in% vars) stop("'vars' must include 'precip' for wet/dry spell analysis", call. = FALSE)

  # date column existence
  if (!"date" %in% names(obs_data[[1]])) stop("'obs_data' must contain a 'date' column", call. = FALSE)
  if (!"date" %in% names(sim_data[[1]][[1]])) stop("'sim_data' must contain a 'date' column", call. = FALSE)

  # variable existence (check only first grid; assumes consistent schema)
  miss <- setdiff(vars, names(obs_data[[1]]))
  if (length(miss) > 0) stop("Variables not found in obs_data: ", paste(miss, collapse = ", "), call. = FALSE)

  if (!is.numeric(wet_q) || length(wet_q) != 1L || !is.finite(wet_q) || wet_q <= 0 || wet_q >= 1) {
    stop("'wet_q' must be between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(extreme_q) || length(extreme_q) != 1L || !is.finite(extreme_q) || extreme_q <= 0 || extreme_q >= 1) {
    stop("'extreme_q' must be between 0 and 1", call. = FALSE)
  }
  if (extreme_q <= wet_q) stop("'extreme_q' must be greater than 'wet_q'", call. = FALSE)

  if (!is.numeric(max_grid) || length(max_grid) != 1L || !is.finite(max_grid) ||
      max_grid < 1 || (max_grid %% 1) != 0) {
    stop("'max_grid' must be a positive integer", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || (seed %% 1) != 0) {
      stop("'seed' must be NULL or a single integer", call. = FALSE)
    }
  }

  invisible(TRUE)
}

# ==============================================================================
# PERIOD STANDARDIZATION (365-day years)
# ==============================================================================

#' Find full (365-day) years after leap-day removal and return longest contiguous block
#' @keywords internal
get_full_year_block <- function(date) {
  if (!inherits(date, "Date")) stop("date must be Date", call. = FALSE)

  yr <- as.integer(format(date, "%Y"))
  tab <- dplyr::tibble(date = date, year = yr) %>%
    dplyr::count(.data$year, name = "n_day") %>%
    dplyr::filter(.data$n_day == 365) %>%
    dplyr::arrange(.data$year)

  if (nrow(tab) == 0) return(list(year = integer(0)))

  years_full <- tab$year
  run_id <- cumsum(c(TRUE, diff(years_full) != 1L))
  runs <- split(years_full, run_id)
  best <- runs[[which.max(vapply(runs, length, integer(1)))]]

  list(year = as.integer(best))
}

#' Pick a random contiguous window of years of length n_year
#' @keywords internal
pick_year_window <- function(year, n_year) {
  year <- as.integer(year)
  if (length(year) < n_year) stop("Not enough full years to pick window.", call. = FALSE)
  if (n_year == length(year)) return(year)

  start_idx <- sample.int(length(year) - n_year + 1L, size = 1L)
  year[start_idx:(start_idx + n_year - 1L)]
}

#' Filter a single grid data.frame to selected years
#' @keywords internal
filter_df_to_years <- function(df, year_keep) {
  if (!("date" %in% names(df))) stop("df must contain 'date'", call. = FALSE)
  yr <- as.integer(format(df$date, "%Y"))
  df[yr %in% year_keep, , drop = FALSE]
}

#' Standardize obs and sim to same full-year length via random windowing
#' @keywords internal
standardize_obs_sim_periods <- function(obs_data, sim_data, n_realizations) {

  # Remove leap days from obs
  obs_leap <- find_leap_day_indices(obs_data[[1]]$date)
  if (!is.null(obs_leap)) obs_data <- lapply(obs_data, function(df) df[-obs_leap, , drop = FALSE])

  # Remove leap days from sim (each realization)
  sim_leap <- find_leap_day_indices(sim_data[[1]][[1]]$date)
  if (!is.null(sim_leap)) {
    sim_data <- lapply(sim_data, function(rlz) lapply(rlz, function(df) df[-sim_leap, , drop = FALSE]))
  }

  obs_years <- get_full_year_block(obs_data[[1]]$date)$year
  sim_years <- get_full_year_block(sim_data[[1]][[1]]$date)$year

  if (length(obs_years) == 0) stop("Observed series has no complete 365-day years after leap-day removal.", call. = FALSE)
  if (length(sim_years) == 0) stop("Simulated series has no complete 365-day years after leap-day removal.", call. = FALSE)

  n_year <- min(length(obs_years), length(sim_years))
  obs_keep <- pick_year_window(obs_years, n_year)
  sim_keep <- pick_year_window(sim_years, n_year)

  obs_data2 <- lapply(obs_data, filter_df_to_years, year_keep = obs_keep)
  sim_data2 <- lapply(seq_len(n_realizations), function(i) lapply(sim_data[[i]], filter_df_to_years, year_keep = sim_keep))

  obs_n <- nrow(obs_data2[[1]])
  sim_n <- nrow(sim_data2[[1]][[1]])
  if (obs_n != sim_n) {
    stop(sprintf("Period standardization failed: obs_n=%d sim_n=%d", obs_n, sim_n), call. = FALSE)
  }

  list(
    obs_data = obs_data2,
    sim_data = sim_data2,
    n_year = n_year,
    obs_year_start = min(obs_keep),
    obs_year_end = max(obs_keep),
    sim_year_start = min(sim_keep),
    sim_year_end = max(sim_keep)
  )
}

# ==============================================================================
# OBS / SIM PROCESSING
# ==============================================================================

#' Process observed weather data
#' @keywords internal
process_observed_data <- function(obs_data, vars, n_grid, wet_q, extreme_q) {

  date <- obs_data[[1]]$date

  obs_key <- dplyr::tibble(
    date = date,
    year = as.integer(format(date, "%Y")),
    mon  = as.integer(format(date, "%m")),
    day  = as.integer(format(date, "%d"))
  )

  obs_long <- lapply(seq_len(n_grid), function(i) {
    df <- obs_data[[i]][, vars, drop = FALSE]
    dplyr::bind_cols(obs_key, df)
  }) %>%
    dplyr::bind_rows(.id = "id") %>%
    dplyr::mutate(id = as.integer(.data$id))

  mc_thresholds <- obs_long %>%
    dplyr::group_by(.data$id, .data$mon) %>%
    dplyr::summarize(
      wet_th = {
        ppos <- .data$precip[is.finite(.data$precip) & .data$precip > 0]
        if (length(ppos) >= 5) stats::quantile(ppos, wet_q, names = FALSE, na.rm = TRUE) else 0
      },
      extreme_th = {
        ppos <- .data$precip[is.finite(.data$precip) & .data$precip > 0]
        if (length(ppos) >= 5) stats::quantile(ppos, extreme_q, names = FALSE, na.rm = TRUE) else 0
      },
      .groups = "drop"
    )

  obs_stats <- compute_timeseries_statistics(
    data = obs_long,
    vars = vars,
    mc_thresholds = mc_thresholds
  )

  list(
    data = obs_long,
    datemat = obs_key,
    mc_thresholds = mc_thresholds,
    stats_season = obs_stats$stats_season %>% dplyr::rename(Observed = .data$value),
    stats_mon_aavg = obs_stats$stats_mon_aavg %>% dplyr::rename(Observed = .data$value),
    stats_annual_aavg = obs_stats$stats_annual_aavg %>% dplyr::rename(Observed = .data$value),
    wetdry = obs_stats$wetdry %>% dplyr::rename(Observed = .data$value),
    cor = obs_stats$cor %>% dplyr::rename(Observed = .data$value),
    cor_cond = obs_stats$cor_cond %>% dplyr::rename(Observed = .data$value)
  )
}

#' Process simulated weather data
#' @keywords internal
process_simulated_data <- function(sim_data, n_realizations, vars, mc_thresholds) {

  sim_key <- dplyr::tibble(
    date = sim_data[[1]][[1]]$date,
    year = as.integer(format(date, "%Y")),
    mon  = as.integer(format(date, "%m")),
    day  = as.integer(format(date, "%d"))
  )

  sim_long <- lapply(seq_len(n_realizations), function(i) {
    sim_data[[i]] %>%
      dplyr::bind_rows(.id = "id") %>%
      dplyr::mutate(id = as.integer(.data$id)) %>%
      dplyr::left_join(sim_key, by = "date")
  })

  sim_stats_list <- lapply(seq_len(n_realizations), function(i) {
    compute_timeseries_statistics(
      data = sim_long[[i]],
      vars = vars,
      mc_thresholds = mc_thresholds
    )
  })

  list(
    stats_season = dplyr::bind_rows(lapply(sim_stats_list, `[[`, "stats_season"), .id = "rlz") %>%
      dplyr::mutate(id = as.numeric(.data$id)) %>%
      dplyr::rename(Simulated = .data$value),

    stats_mon_aavg = dplyr::bind_rows(lapply(sim_stats_list, `[[`, "stats_mon_aavg"), .id = "rlz") %>%
      dplyr::rename(Simulated = .data$value),

    stats_annual_aavg = dplyr::bind_rows(lapply(sim_stats_list, `[[`, "stats_annual_aavg"), .id = "rlz") %>%
      dplyr::mutate(year = .data$year - min(.data$year) + 1) %>%
      dplyr::rename(Simulated = .data$value),

    cor = dplyr::bind_rows(lapply(sim_stats_list, `[[`, "cor"), .id = "rlz") %>%
      dplyr::rename(Simulated = .data$value),

    wetdry = dplyr::bind_rows(lapply(sim_stats_list, `[[`, "wetdry"), .id = "rlz") %>%
      dplyr::mutate(id = as.numeric(.data$id)) %>%
      dplyr::rename(Simulated = .data$value),

    cor_cond = dplyr::bind_rows(lapply(sim_stats_list, `[[`, "cor_cond"), .id = "rlz") %>%
      dplyr::rename(Simulated = .data$value)
  )
}

# ==============================================================================
# STATISTICS CORE (aligned to precip)
# ==============================================================================

#' Compute fit metrics for each realization
#'
#' Computes Mean Absolute Error (MAE) for all key metrics.
#' MAE measures the average magnitude of errors between simulated and observed values.
#' Lower MAE indicates better fit.
#'
#' @param obs_results Observed data results from process_observed_data
#' @param sim_results Simulated data results from process_simulated_data
#' @param variables Character vector of variable names
#'
#' @keywords internal
compute_realization_fit_metrics <- function(obs_results, sim_results, vars) {

  metrics_list <- list()

  # --- 1. MAE for means (by variable, across all grid-months) ---
  mae_mean <- sim_results$stats_season %>%
    dplyr::filter(.data$stat == "mean") %>%
    dplyr::left_join(
      obs_results$stats_season %>% dplyr::select(.data$id, .data$mon, .data$variable, Observed),
      by = c("id", "mon", "variable")
    ) %>%
    dplyr::group_by(.data$rlz, .data$variable) %>%
    dplyr::summarize(
      mae_mean = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$mae_mean, names_prefix = "mae_mean_")

  # --- 2. MAE for SDs ---
  mae_sd <- sim_results$stats_season %>%
    dplyr::filter(.data$stat == "sd") %>%
    dplyr::left_join(
      obs_results$stats_season %>% dplyr::select(.data$id, .data$mon, .data$variable, Observed),
      by = c("id", "mon", "variable")
    ) %>%
    dplyr::group_by(.data$rlz, .data$variable) %>%
    dplyr::summarize(
      mae_sd = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$mae_sd, names_prefix = "mae_sd_")

  # --- 3. MAE for wet/dry day counts ---
  mae_wetdry_days <- sim_results$wetdry %>%
    dplyr::filter(.data$type == "days") %>%
    dplyr::left_join(
      obs_results$wetdry %>% dplyr::select(.data$id, .data$mon, .data$stat, Observed),
      by = c("id", "mon", "stat")
    ) %>%
    dplyr::group_by(.data$rlz, .data$stat) %>%
    dplyr::summarize(
      mae = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$stat, values_from = .data$mae, names_prefix = "mae_days_")

  # --- 4. MAE for spell lengths ---
  mae_spells <- sim_results$wetdry %>%
    dplyr::filter(.data$type == "spells") %>%
    dplyr::left_join(
      obs_results$wetdry %>% dplyr::select(.data$id, .data$mon, .data$stat, Observed),
      by = c("id", "mon", "stat")
    ) %>%
    dplyr::group_by(.data$rlz, .data$stat) %>%
    dplyr::summarize(
      mae = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$stat, values_from = .data$mae, names_prefix = "mae_spell_")

  # --- 5. MAE for cross-grid correlations ---
  mae_cor_crossgrid <- sim_results$cor %>%
    dplyr::filter(.data$variable1 == .data$variable2, .data$id1 != .data$id2) %>%
    dplyr::left_join(
      obs_results$cor %>%
        dplyr::filter(.data$variable1 == .data$variable2, .data$id1 != .data$id2) %>%
        dplyr::select(.data$id1, .data$variable1, .data$id2, .data$variable2, Observed),
      by = c("id1", "variable1", "id2", "variable2")
    ) %>%
    dplyr::group_by(.data$rlz) %>%
    dplyr::summarize(
      mae_cor_crossgrid = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    )

  # --- 6. MAE for inter-variable correlations ---
  mae_cor_intervariable <- sim_results$cor %>%
    dplyr::filter(.data$id1 == .data$id2, .data$variable1 != .data$variable2) %>%
    dplyr::left_join(
      obs_results$cor %>%
        dplyr::filter(.data$id1 == .data$id2, .data$variable1 != .data$variable2) %>%
        dplyr::select(.data$id1, .data$variable1, .data$id2, .data$variable2, Observed),
      by = c("id1", "variable1", "id2", "variable2")
    ) %>%
    dplyr::group_by(.data$rlz) %>%
    dplyr::summarize(
      mae_cor_intervariable = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    )

  # --- 7. Merge all metrics ---
  summary_df <- mae_mean %>%
    dplyr::left_join(mae_sd, by = "rlz") %>%
    dplyr::left_join(mae_wetdry_days, by = "rlz") %>%
    dplyr::left_join(mae_spells, by = "rlz") %>%
    dplyr::left_join(mae_cor_crossgrid, by = "rlz") %>%
    dplyr::left_join(mae_cor_intervariable, by = "rlz")

  # --- 8. Compute overall score (normalized, lower = better) ---
  # Normalize all MAE metrics to 0-1 scale for fair weighting when combining
  metric_cols <- setdiff(names(summary_df), "rlz")

  summary_df <- summary_df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(metric_cols),
        ~ {
          rng <- range(.x, na.rm = TRUE)
          den <- rng[2] - rng[1]
          if (!is.finite(den) || den == 0) {
            # constant / all-NA column -> neutral 0 (won't affect means with na.rm=TRUE)
            return(rep(0, length(.x)))
          }
          (.x - rng[1]) / den
        },
        .names = "norm_{.col}"
      )
    )

  norm_cols <- grep("^norm_", names(summary_df), value = TRUE)

  summary_df <- summary_df %>%
    dplyr::mutate(
      overall_score = rowMeans(dplyr::across(dplyr::all_of(norm_cols)), na.rm = TRUE)
    ) %>%
    dplyr::select(-dplyr::all_of(norm_cols)) %>%
    dplyr::arrange(.data$overall_score) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::select(.data$rlz, .data$rank, .data$overall_score, dplyr::everything())

  summary_df
}

#' Compute time series statistics used by the assessment
#' @keywords internal
compute_timeseries_statistics <- function(data, vars, mc_thresholds) {

  stat_funs <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)
  year_n <- length(unique(format(data$date, "%Y")))

  stats_season <- compute_grouped_statistics(
    df = data,
    vars = vars,
    group_var = c("id", "mon"),
    stat_funs = stat_funs
  )

  stats_mon_aavg <- compute_grouped_statistics(
    df = data,
    vars = vars,
    group_var = c("year", "mon"),
    stat_funs = stat_funs
  )

  stats_annual_aavg <- compute_grouped_statistics(
    df = data,
    vars = vars,
    group_var = "year",
    stat_funs = stat_funs
  ) %>%
    dplyr::mutate(year = .data$year - min(.data$year) + 1)

  wetdry <- data %>%
    dplyr::left_join(mc_thresholds, by = c("id", "mon")) %>%
    dplyr::group_by(.data$id, .data$mon) %>%
    dplyr::summarize(
      Wet_days = {
        th <- .data$wet_th[1]
        if (!is.finite(th) || th <= 0) sum(.data$precip > 0, na.rm = TRUE) / year_n
        else sum(.data$precip > th, na.rm = TRUE) / year_n
      },
      Dry_days = {
        th <- .data$wet_th[1]
        if (!is.finite(th) || th <= 0) sum(.data$precip <= 0, na.rm = TRUE) / year_n
        else sum(.data$precip <= th, na.rm = TRUE) / year_n
      },
      Dry_spells = {
        th <- .data$wet_th[1]
        thr <- if (!is.finite(th) || th <= 0) 0 else th
        lens <- compute_spell_lengths(.data$precip, threshold = thr, below = TRUE)
        if (length(lens) == 0) 0 else mean(lens)
      },
      Wet_spells = {
        th <- .data$wet_th[1]
        thr <- if (!is.finite(th) || th <= 0) 0 else th
        lens <- compute_spell_lengths(.data$precip, threshold = thr, below = FALSE)
        if (length(lens) == 0) 0 else mean(lens)
      },
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(.data$Wet_days, .data$Dry_days, .data$Dry_spells, .data$Wet_spells),
      names_to = "stat_full",
      values_to = "value"
    ) %>%
    tidyr::separate(.data$stat_full, into = c("stat", "type"), sep = "_") %>%
    dplyr::mutate(variable = "precip", .after = .data$mon)

  cor <- compute_correlation_matrix_anom(data, vars)

  cor_cond <- compute_conditional_precip_cor(
    data = data,
    vars = vars,
    mc_thresholds = mc_thresholds,
    wet_def = "monthly_quantile",
    use_anom = TRUE,
    use_log_precip_on_wet = TRUE,
    method = "pearson",
    min_pairs = 50
  )

  list(
    stats_season = stats_season,
    stats_mon_aavg = stats_mon_aavg,
    stats_annual_aavg = stats_annual_aavg,
    wetdry = wetdry,
    cor = cor,
    cor_cond = cor_cond
  )
}

#' Compute grouped statistics (long-form output)
#' @keywords internal
compute_grouped_statistics <- function(df, vars, group_var, stat_funs) {

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
    dplyr::summarize(
      dplyr::across(dplyr::all_of(vars), stat_funs, .names = "{.col}:{.fn}"),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(group_var),
      names_to = "variable",
      values_to = "value"
    ) %>%
    tidyr::separate(.data$variable, into = c("variable", "stat"), sep = ":") %>%
    dplyr::mutate(
      stat = factor(.data$stat, levels = names(stat_funs)),
      variable = as.character(.data$variable)
    )
}

#' Compute correlation matrix on monthly anomalies (grid:variable pairs)
#' @keywords internal
compute_correlation_matrix_anom <- function(data, vars) {

  dat2 <- data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::group_by(.data$id, .data$mon, .data$variable) %>%
    dplyr::mutate(value = value - mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    tidyr::unite("id_variable", .data$id, .data$variable, sep = ":") %>%
    dplyr::select(.data$date, .data$id_variable, .data$value) %>%
    tidyr::pivot_wider(names_from = .data$id_variable, values_from = .data$value) %>%
    dplyr::arrange(.data$date)

  mat <- as.matrix(dat2 %>% dplyr::select(-.data$date))
  cmat <- stats::cor(mat, use = "pairwise.complete.obs", method = "pearson")
  tri <- upper.tri(cmat, diag = FALSE)

  dplyr::tibble(
    id_variable1 = colnames(cmat)[row(cmat)[tri]],
    id_variable2 = colnames(cmat)[col(cmat)[tri]],
    value = cmat[tri]
  ) %>%
    tidyr::separate(.data$id_variable1, c("id1", "variable1"), sep = ":") %>%
    tidyr::separate(.data$id_variable2, c("id2", "variable2"), sep = ":")
}

#' Conditional precip–X correlations within each grid
#' @keywords internal
compute_conditional_precip_cor <- function(data,
                                         vars,
                                         mc_thresholds = NULL,
                                         wet_def = c("gt0", "monthly_quantile"),
                                         use_anom = TRUE,
                                         use_log_precip_on_wet = TRUE,
                                         method = "pearson",
                                         min_pairs = 50) {

  wet_def <- match.arg(wet_def)

  if (!("precip" %in% vars)) stop("vars must include 'precip'", call. = FALSE)
  targets <- setdiff(vars, "precip")
  if (length(targets) == 0) {
    return(dplyr::tibble(
      id1 = integer(), variable1 = character(),
      id2 = integer(), variable2 = character(),
      regime = character(), transform = character(),
      value = numeric(), n = integer()
    ))
  }

  dat <- data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("precip", targets)),
      names_to = "variable",
      values_to = "value"
    )

  if (use_anom) {
    dat <- dat %>%
      dplyr::group_by(.data$id, .data$mon, .data$variable) %>%
      dplyr::mutate(value = value - mean(value, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }

  dat_wide <- dat %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$value)

  if (wet_def == "monthly_quantile") {
    if (is.null(mc_thresholds)) stop("mc_thresholds required when wet_def='monthly_quantile'", call. = FALSE)
    dat_wide <- dat_wide %>% dplyr::left_join(mc_thresholds, by = c("id", "mon"))
  }

  .safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if (n < min_pairs) return(c(value = NA_real_, n = n))
    c(value = stats::cor(x[ok], y[ok], method = method), n = n)
  }

  by_id <- split(dat_wide, dat_wide$id)

  out <- lapply(by_id, function(df_id) {

    wet_flag <- if (wet_def == "monthly_quantile") {
      df_id$precip >= df_id$wet_th
    } else {
      df_id$precip > 0
    }
    dry_flag <- !wet_flag & !is.na(df_id$precip)

    lapply(targets, function(v) {

      res_all <- .safe_cor(df_id$precip, df_id[[v]])

      if (use_log_precip_on_wet) {
        p_wet <- log1p(pmax(df_id$precip, 0))
        res_wet <- .safe_cor(p_wet[wet_flag], df_id[[v]][wet_flag])
        tr_wet <- "log1p_precip"
      } else {
        res_wet <- .safe_cor(df_id$precip[wet_flag], df_id[[v]][wet_flag])
        tr_wet <- "precip"
      }

      res_dry <- .safe_cor(df_id$precip[dry_flag], df_id[[v]][dry_flag])

      dplyr::tibble(
        id1 = df_id$id[1],
        variable1 = "precip",
        id2 = df_id$id[1],
        variable2 = v,
        regime = c("all", "wet", "dry"),
        transform = c("precip", tr_wet, "precip"),
        value = c(res_all["value"], res_wet["value"], res_dry["value"]),
        n = c(as.integer(res_all["n"]), as.integer(res_wet["n"]), as.integer(res_dry["n"]))
      )
    }) %>% dplyr::bind_rows()

  }) %>% dplyr::bind_rows()

  out
}


#' Filter correlations by type (cross-grid or inter-variable)
#'
#' Internal helper to subset correlation data into the intended pairing type:
#' - cross-grid correlations: same variable, different grid ids
#' - inter-variable correlations: different variables, same grid id
#'
#' @param cor_data data.frame. Correlation data with columns id1, id2, variable1, variable2.
#' @param same_var Logical. Keep rows where variable1 == variable2 (TRUE) or != (FALSE).
#' @param same_id Logical. Keep rows where id1 == id2 (TRUE) or != (FALSE).
#' @param allowed_pairs Character vector of order-invariant pair keys.
#'   For pair_type="id": "min(id1,id2):max(id1,id2)"
#'   For pair_type="variable": "min(var1,var2):max(var1,var2)"
#' @param pair_type Character. Either "id" or "variable".
#'
#' @return Filtered data.frame.
#' @keywords internal
filter_correlations <- function(cor_data, same_var, same_id, allowed_pairs, pair_type) {

  if (!is.data.frame(cor_data)) stop("cor_data must be a data.frame", call. = FALSE)
  req <- c("id1", "id2", "variable1", "variable2")
  miss <- setdiff(req, names(cor_data))
  if (length(miss) > 0) {
    stop("cor_data is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!pair_type %in% c("id", "variable")) stop("pair_type must be 'id' or 'variable'", call. = FALSE)

  out <- cor_data %>%
    dplyr::filter(
      (.data$variable1 == .data$variable2) == same_var,
      (.data$id1 == .data$id2) == same_id
    )

  if (length(out) == 0) return(out)

  if (pair_type == "id") {
    out <- out %>%
      dplyr::mutate(
        pair = paste(pmin(.data$id1, .data$id2), pmax(.data$id1, .data$id2), sep = ":")
      ) %>%
      dplyr::filter(.data$pair %in% allowed_pairs) %>%
      dplyr::select(-.data$pair)
  } else {
    out <- out %>%
      dplyr::mutate(
        pair = paste(
          pmin(.data$variable1, .data$variable2),
          pmax(.data$variable1, .data$variable2),
          sep = ":"
        ),
        variable = .data$pair
      ) %>%
      dplyr::filter(.data$pair %in% allowed_pairs) %>%
      dplyr::select(-.data$pair)
  }

  out
}


#' Display fit summary table for all realizations
#' @keywords internal
display_fit_summary_table <- function(fit_summary, variables) {

  if (is.null(fit_summary) || nrow(fit_summary) == 0) {
    cat("\n")
    cat("WARNING: No fit metrics computed\n")
    cat("\n")
    return(invisible(NULL))
  }

  # Select key metrics for display (good coverage of different features)
  key_metrics <- character(0)

  # 1. Mean accuracy - select first variable if precip exists, otherwise first variable
  if ("mae_mean_precip" %in% names(fit_summary)) {
    key_metrics <- c(key_metrics, "mae_mean_precip")
  } else {
    mean_cols <- grep("^mae_mean_", names(fit_summary), value = TRUE)
    if (length(mean_cols) > 0) key_metrics <- c(key_metrics, mean_cols[1])
  }

  # 2. Variability accuracy - select first variable
  sd_cols <- grep("^mae_sd_", names(fit_summary), value = TRUE)
  if (length(sd_cols) > 0) key_metrics <- c(key_metrics, sd_cols[1])

  # 3. Wet/dry day counts
  if ("mae_days_Wet" %in% names(fit_summary)) {
    key_metrics <- c(key_metrics, "mae_days_Wet")
  }

  # 4. Spell lengths
  if ("mae_spell_Wet" %in% names(fit_summary)) {
    key_metrics <- c(key_metrics, "mae_spell_Wet")
  }

  # 5. Cross-grid correlation
  if ("mae_cor_crossgrid" %in% names(fit_summary)) {
    key_metrics <- c(key_metrics, "mae_cor_crossgrid")
  }

  # 6. Inter-variable correlation
  if ("mae_cor_intervariable" %in% names(fit_summary)) {
    key_metrics <- c(key_metrics, "mae_cor_intervariable")
  }

  # 7. Overall score (always include)
  key_metrics <- c(key_metrics, "overall_score")

  # Filter to available metrics
  key_metrics <- key_metrics[key_metrics %in% names(fit_summary)]

  # Prepare display data
  display_data <- fit_summary[, c("rlz", "rank", key_metrics), drop = FALSE]

  # Create nice column names for display
  col_names <- names(display_data)
  display_names <- col_names
  display_names <- gsub("^mae_mean_", "Mean.", display_names)
  display_names <- gsub("^mae_sd_", "SD.", display_names)
  display_names <- gsub("^mae_days_", "Days.", display_names)
  display_names <- gsub("^mae_spell_", "Spell.", display_names)
  display_names <- gsub("mae_cor_crossgrid", "Cor.Cross", display_names)
  display_names <- gsub("mae_cor_intervariable", "Cor.Inter", display_names)
  display_names <- gsub("overall_score", "Score", display_names)
  display_names <- gsub("rlz", "Rlz", display_names)
  display_names <- gsub("rank", "Rank", display_names)

  # Calculate column widths
  col_widths <- pmax(
    nchar(display_names),
    apply(display_data, 2, function(x) max(nchar(format(x, digits = 4, nsmall = 4))))
  ) + 2

  # Ensure minimum widths
  col_widths <- pmax(col_widths, 6)

  # Print header
  cat("\n")
  cat(strrep("=", sum(col_widths) + length(col_widths) - 1), "\n")
  cat(" FIT ASSESSMENT SUMMARY - ALL REALIZATIONS\n")
  cat(strrep("=", sum(col_widths) + length(col_widths) - 1), "\n")
  cat("\n")

  # Print metric descriptions
  cat(" Metrics shown (Mean Absolute Error, lower is better):\n")
  cat("  MAE measures average magnitude of errors between simulated and observed\n")
  cat("\n")
  for (i in seq_along(key_metrics)) {
    if (key_metrics[i] == "overall_score") {
      cat("  - Score       : Normalized overall score across all metrics\n")
    } else if (grepl("^mae_mean_", key_metrics[i])) {
      var <- sub("^mae_mean_", "", key_metrics[i])
      cat("  - Mean.", var, "   : MAE of monthly means\n", sep = "")
    } else if (grepl("^mae_sd_", key_metrics[i])) {
      var <- sub("^mae_sd_", "", key_metrics[i])
      cat("  - SD.", var, "     : MAE of monthly standard deviations\n", sep = "")
    } else if (grepl("^mae_days_", key_metrics[i])) {
      cat("  - Days.", sub("^mae_days_", "", key_metrics[i]), "   : MAE of wet/dry day counts\n", sep = "")
    } else if (grepl("^mae_spell_", key_metrics[i])) {
      cat("  - Spell.", sub("^mae_spell_", "", key_metrics[i]), "  : MAE of wet/dry spell lengths (days)\n", sep = "")
    } else if (key_metrics[i] == "mae_cor_crossgrid") {
      cat("  - Cor.Cross   : MAE of cross-grid correlations\n")
    } else if (key_metrics[i] == "mae_cor_intervariable") {
      cat("  - Cor.Inter   : MAE of inter-variable correlations\n")
    }
  }
  cat("\n")

  # Print column headers
  header_line <- paste(
    mapply(function(name, width) {
      format(name, width = width, justify = "right")
    }, display_names, col_widths),
    collapse = " "
  )
  cat(header_line, "\n")
  cat(strrep("-", sum(col_widths) + length(col_widths) - 1), "\n")

  # Print data rows
  for (i in seq_len(nrow(display_data))) {
    row_values <- as.character(display_data[i, ])

    # Format numeric values
    for (j in seq_along(row_values)) {
      if (col_names[j] %in% c("rlz", "rank")) {
        row_values[j] <- format(as.integer(display_data[i, j]), width = col_widths[j], justify = "right")
      } else {
        row_values[j] <- format(
          round(as.numeric(display_data[i, j]), 4),
          width = col_widths[j],
          justify = "right",
          nsmall = 4
        )
      }
    }

    cat(paste(row_values, collapse = " "), "\n")
  }

  cat(strrep("=", sum(col_widths) + length(col_widths) - 1), "\n")

  # Print summary statistics
  best_rlz <- fit_summary$rlz[1]
  worst_rlz <- fit_summary$rlz[nrow(fit_summary)]
  median_score <- median(fit_summary$overall_score, na.rm = TRUE)

  cat("\n")
  cat(" Summary:\n")
  cat("  - Best realization  : ", best_rlz, " (score = ",
      format(round(fit_summary$overall_score[1], 4), nsmall = 4), ")\n", sep = "")
  cat("  - Worst realization : ", worst_rlz, " (score = ",
      format(round(fit_summary$overall_score[nrow(fit_summary)], 4), nsmall = 4), ")\n", sep = "")
  cat("  - Median score      : ", format(round(median_score, 4), nsmall = 4), "\n", sep = "")
  cat("\n")

  invisible(fit_summary)
}


# ==============================================================================
# PLOT DATA PREP
# ==============================================================================

#' Prepare data for diagnostic plotting (obs vs sim)
#' @keywords internal
prepare_plot_data <- function(obs_results, sim_results, vars) {

  stat_funs <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)

  daily_stats_season <- sim_results$stats_season %>%
    dplyr::left_join(obs_results$stats_season, by = c("id", "mon", "variable", "stat")) %>%
    dplyr::mutate(
      variable = factor(.data$variable),
      stat = factor(.data$stat, levels = names(stat_funs))
    )

  stats_allcor <- sim_results$cor %>%
    dplyr::left_join(obs_results$cor, by = c("id1", "variable1", "id2", "variable2"))

  stats_allcor_cond <- sim_results$cor_cond %>%
    dplyr::left_join(obs_results$cor_cond, by = c("id1", "variable1", "id2", "variable2", "regime", "transform"))

  ids <- sort(unique(obs_results$data$id))
  id_pairs_allowed <- if (length(ids) >= 2) {
    apply(utils::combn(ids, 2), 2, function(x) paste(x[1], x[2], sep = ":"))
  } else character(0)

  vars_use <- sort(unique(as.character(vars)))
  var_pairs_allowed <- if (length(vars_use) >= 2) {
    apply(utils::combn(vars_use, 2), 2, function(x) paste(x[1], x[2], sep = ":"))
  } else character(0)

  stats_crosscor <- filter_correlations(
    cor_data = stats_allcor,
    same_var = TRUE,
    same_id = FALSE,
    allowed_pairs = id_pairs_allowed,
    pair_type = "id"
  )

  stats_intercor <- filter_correlations(
    cor_data = stats_allcor,
    same_var = FALSE,
    same_id = TRUE,
    allowed_pairs = var_pairs_allowed,
    pair_type = "variable"
  )

  stats_wetdry <- sim_results$wetdry %>%
    dplyr::left_join(obs_results$wetdry, by = c("id", "mon", "variable", "stat", "type")) %>%
    dplyr::mutate(stat = factor(.data$stat, levels = c("Dry", "Wet")))

  list(
    daily_stats_season = daily_stats_season,
    stats_mon_aavg_sim = sim_results$stats_mon_aavg,
    stats_mon_aavg_obs = obs_results$stats_mon_aavg,
    stats_annual_aavg_sim = sim_results$stats_annual_aavg,
    stats_annual_aavg_obs = obs_results$stats_annual_aavg,
    stats_crosscor = stats_crosscor,
    stats_intercor = stats_intercor,
    stats_wetdry = stats_wetdry,
    stats_precip_cor_cond = stats_allcor_cond
  )
}
