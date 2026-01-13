#' Evaluate Stochastic Weather Generator Performance
#'
#' Comprehensive diagnostic evaluation comparing synthetic weather simulations
#' against historical observations across multiple grid cells. Computes summary
#' statistics, wet/dry day counts, spell lengths, and inter-site correlations,
#' and generates diagnostic plots to assess stochastic weather generator performance.
#'
#' @param daily.sim List of simulated weather realizations. Each element should be
#'   a list of data frames (one per grid cell), containing daily values and a `date` column.
#' @param daily.obs List of observed weather data frames (one per grid cell). Each
#'   should contain a `date` column and the variables specified in `variables`.
#' @param variables Character vector of variable names to evaluate (e.g., `c("precip", "temp")`).
#'   Must include "precip" for wet/dry spell analysis.
#' @param variable.labels Optional character vector of variable labels for plots.
#'   Defaults to `variables` if `NULL`.
#' @param realization.num Integer. Number of synthetic realizations in `daily.sim`.
#' @param wet.quantile Numeric between 0 and 1. Quantile threshold for wet days
#'   (default = 0.2).
#' @param extreme.quantile Numeric between 0 and 1. Quantile threshold for extremely
#'   wet days (default = 0.8).
#' @param output.path Character. Directory path to save generated plots. If `NULL`,
#'   plots are not saved to disk.
#' @param save.plots Logical. Whether to save plots to `output.path` (default = `TRUE`).
#' @param show.title Logical. Whether to display titles in plots (default = `TRUE`).
#' @param max.grids Integer. Maximum number of grid cells to evaluate (default = 25).
#'   If more grids are provided, a random subsample is used to control memory.
#' @param seed Optional integer. Random seed for reproducible grid subsampling and
#'   year window selection. If NULL, results will vary between runs.
#'
#' @return A named list of `ggplot2` plot objects with class "weather_assessment".
#'   The returned object also contains fit summary metrics stored as attributes.
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom e1071 skewness
#' @importFrom stats cor sd quantile
#' @importFrom utils combn
#' @export
evaluate_weather_generator <- function(
    daily.sim = NULL,
    daily.obs = NULL,
    variables = NULL,
    variable.labels = NULL,
    realization.num = NULL,
    wet.quantile = 0.2,
    extreme.quantile = 0.8,
    output.path = NULL,
    save.plots = TRUE,
    show.title = TRUE,
    max.grids = 25,
    seed = NULL
) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  validate_inputs(
    daily.sim = daily.sim,
    daily.obs = daily.obs,
    variables = variables,
    realization.num = realization.num,
    wet.quantile = wet.quantile,
    extreme.quantile = extreme.quantile
  )

  if (!is.numeric(max.grids) ||
      length(max.grids) != 1L ||
      is.na(max.grids) ||
      !is.finite(max.grids) ||
      max.grids < 1 ||
      max.grids != as.integer(max.grids)) {
    stop("'max.grids' must be a positive integer", call. = FALSE)
  }

  # ============================================================================
  # RNG STATE MANAGEMENT
  # ============================================================================

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or a single finite number", call. = FALSE)
    }
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }
    on.exit({ if (has_seed) .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed)
  }

  # ============================================================================
  # SETUP
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(
      paste0(
        "[Assessment] Start | grids = {length(daily.obs)} | ",
        "realizations = {realization.num} | ")
    )
  }

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(
      paste0(
        "[Assessment] variables = {paste(variables, collapse = ',')} | ")
    )
  }

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(
      paste0("[Assessment] Parameters: wet.q = { wet.quantile} | extreme.q = { extreme.quantile}")
    )
  }


  if (is.null(variable.labels)) variable.labels <- variables

  if (!is.null(output.path)) {
    if (!dir.exists(output.path)) {
      dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    save.plots <- FALSE
  }

  options(dplyr.summarise.inform = FALSE, tidyverse.quiet = TRUE)

  plot.config <- list(
    subtitle = "Value range and median from all simulations shown against observed",
    alpha = 0.4,
    colors = stats::setNames(c("blue3", "gray40"), c("Observed", "Simulated")),
    theme = ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14),
        plot.subtitle = ggplot2::element_text(size = 10)
      )
  )

  plots <- list()

  # ============================================================================
  # GRID SUBSAMPLING TO CONTROL MEMORY USE
  # ============================================================================

  n_grids <- length(daily.obs)
  n_grids_org <- n_grids

  if (n_grids > max.grids) {

    sel.grids <- sort(sample(seq_len(n_grids), max.grids))

    daily.obs <- daily.obs[sel.grids]
    daily.sim <- lapply(daily.sim, function(rlz) rlz[sel.grids])

    n_grids <- length(daily.obs)

    if (requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn(
        "[Assessment] Grid count reduced from {n_grids_org} to {n_grids} for memory control"
      )
    }
  }


  # ============================================================================
  # STANDARDIZE PERIODS (FULL YEARS + MATCH LENGTH VIA RANDOM WINDOW)
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Standardizing obs/sim periods to full years and equal length")
  }

  std <- standardize_obs_sim_periods(
    daily.obs = daily.obs,
    daily.sim = daily.sim,
    realization.num = realization.num,
    variables = variables
  )

  daily.obs <- std$daily.obs
  daily.sim <- std$daily.sim

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(
      paste0(
        "[Assessment] Standardized period | ",
        "Obs = ", std$obs_year_start, "-", std$obs_year_end, " | ",
        "Sim = ", std$sim_year_start, "-", std$sim_year_end
      )
    )
  }

  # ============================================================================
  # PROCESS OBSERVED DATA
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Processing observed data")
  }

  obs.results <- process_observed_data(
    daily.obs = daily.obs,
    variables = variables,
    n_grids = n_grids,
    wet.quantile = wet.quantile,
    extreme.quantile = extreme.quantile
  )

  # ============================================================================
  # PROCESS SIMULATED DATA
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Processing simulated data ({realization.num} realizations)")
  }

  sim.results <- process_simulated_data(
    daily.sim = daily.sim,
    realization.num = realization.num,
    variables = variables,
    mc.thresholds = obs.results$mc.thresholds
  )

  # ============================================================================
  # MERGE AND PREPARE PLOT DATA
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Preparing diagnostic data for plotting")
  }

  plot.data <- prepare_plot_data(
    obs.results = obs.results,
    sim.results = sim.results,
    variables = variables,
    n_grids = n_grids
  )

  # ============================================================================
  # GENERATE DIAGNOSTIC PLOTS
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Generating diagnostic plots")
  }

  plots <- create_all_diagnostic_plots(
    plot.data = plot.data,
    plot.config = plot.config,
    variables = variables,
    variable.labels = variable.labels,
    show.title = show.title,
    save.plots = save.plots,
    output.path = output.path
  )

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Generated {length(plots)} diagnostic plots")
    if (save.plots) {
      logger::log_info("[Assessment] Plots saved to: {output.path}")
    }
  }

  # ============================================================================
  # COMPUTE FIT METRICS SUMMARY TABLE
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Computing fit metrics for all realizations")
  }

  fit_summary <- compute_realization_fit_metrics(
    obs.results = obs.results,
    sim.results = sim.results,
    variables = variables
  )

  # ============================================================================
  # DISPLAY FIT SUMMARY TABLE
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Displaying fit assessment summary")
  }

  display_fit_summary_table(fit_summary, variables)

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Assessment completed successfully")
  }

  structure(
    plots,
    class = c("weather_assessment", "list"),
    fit_summary = fit_summary,
    metadata = list(
      n_grids = n_grids,
      realization.num = realization.num,
      variables = variables,
      assessment.date = Sys.Date()
    )
  )
}

# ==============================================================================
# FIT SUMMARY DISPLAY
# ==============================================================================

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
# HELPERS
# ==============================================================================

#' Find full (365-day) years after leap-day removal and return longest contiguous block
#' @keywords internal
get_full_year_block <- function(dates) {

  if (!inherits(dates, "Date")) stop("dates must be Date", call. = FALSE)

  yrs <- as.integer(format(dates, "%Y"))
  tab <- dplyr::tibble(date = dates, year = yrs) %>%
    dplyr::count(.data$year, name = "n_days") %>%
    dplyr::filter(.data$n_days == 365) %>%
    dplyr::arrange(.data$year)

  if (nrow(tab) == 0) {
    return(list(years = integer(0)))
  }

  years_full <- tab$year
  # Build contiguous runs
  run_id <- cumsum(c(TRUE, diff(years_full) != 1L))
  runs <- split(years_full, run_id)

  # Take the longest run; if tie, take the first
  lens <- vapply(runs, length, integer(1))
  best <- runs[[which.max(lens)]]

  list(years = as.integer(best))
}

#' Pick a random contiguous window of years of length n_years from a year vector
#' @keywords internal
pick_year_window <- function(years, n_years) {

  years <- as.integer(years)
  if (length(years) < n_years) stop("Not enough full years to pick window.", call. = FALSE)
  if (n_years == length(years)) return(years)

  start_idx <- sample.int(length(years) - n_years + 1L, size = 1L)
  years[start_idx:(start_idx + n_years - 1L)]
}

#' Filter a single grid df to selected years (keeps original columns)
#' @keywords internal
filter_df_to_years <- function(df, years_keep) {

  if (!("date" %in% names(df))) stop("df must contain 'date'", call. = FALSE)

  yrs <- as.integer(format(df$date, "%Y"))
  df[yrs %in% years_keep, , drop = FALSE]
}

#' Standardize obs and sim to same full-year length via random windowing
#' Applies the chosen window consistently across all grids and all realizations.
#' @keywords internal
standardize_obs_sim_periods <- function(daily.obs, daily.sim, realization.num, variables) {

  # 1) Remove leap days consistently (based on each series' date vector)
  #    Assumption: all grids share the same date vector within obs and within each realization.
  obs_dates <- daily.obs[[1]]$date
  sim_dates <- daily.sim[[1]][[1]]$date

  obs_leap <- find_leap_days(obs_dates)
  sim_leap <- find_leap_days(sim_dates)

  if (!is.null(obs_leap)) {
    daily.obs <- lapply(daily.obs, function(df) df[-obs_leap, , drop = FALSE])
    obs_dates <- daily.obs[[1]]$date
  }

  if (!is.null(sim_leap)) {
    daily.sim <- lapply(daily.sim, function(rlz) {
      lapply(rlz, function(df) df[-sim_leap, , drop = FALSE])
    })
    sim_dates <- daily.sim[[1]][[1]]$date
  }

  # 2) Identify longest contiguous full-year blocks
  obs_block <- get_full_year_block(obs_dates)$years
  sim_block <- get_full_year_block(sim_dates)$years

  if (length(obs_block) == 0) stop("Observed series has no complete 365-day years after leap-day removal.", call. = FALSE)
  if (length(sim_block) == 0) stop("Simulated series has no complete 365-day years after leap-day removal.", call. = FALSE)

  n_years <- min(length(obs_block), length(sim_block))

  # 3) Pick windows (randomly within each block if longer than needed)
  obs_years_keep <- pick_year_window(obs_block, n_years)
  sim_years_keep <- pick_year_window(sim_block, n_years)

  # 4) Apply year filtering
  daily.obs2 <- lapply(daily.obs, filter_df_to_years, years_keep = obs_years_keep)

  daily.sim2 <- lapply(seq_len(realization.num), function(i) {
    lapply(daily.sim[[i]], filter_df_to_years, years_keep = sim_years_keep)
  })

  # Basic integrity check: after filtering, all grids should have same nrow within obs and sim
  obs_n <- nrow(daily.obs2[[1]])
  sim_n <- nrow(daily.sim2[[1]][[1]])

  if (obs_n != sim_n) {
    stop(
      "Period standardization failed: obs and sim window lengths differ after filtering. ",
      "obs_n=", obs_n, " sim_n=", sim_n,
      call. = FALSE
    )
  }

  list(
    daily.obs = daily.obs2,
    daily.sim = daily.sim2,
    n_years = n_years,
    obs_year_start = min(obs_years_keep),
    obs_year_end = max(obs_years_keep),
    sim_year_start = min(sim_years_keep),
    sim_year_end = max(sim_years_keep)
  )
}

#' Validate inputs for weather assessment
#' @keywords internal
validate_inputs <- function(daily.sim, daily.obs, variables, realization.num,
                            wet.quantile, extreme.quantile) {

  if (is.null(daily.sim)) stop("'daily.sim' must not be NULL")
  if (is.null(daily.obs)) stop("'daily.obs' must not be NULL")
  if (is.null(variables)) stop("'variables' must not be NULL")
  if (is.null(realization.num)) stop("'realization.num' must not be NULL")

  if (!is.list(daily.sim)) stop("'daily.sim' must be a list")
  if (!is.list(daily.obs)) stop("'daily.obs' must be a list")
  if (!is.character(variables)) stop("'variables' must be a character vector")
  if (!is.numeric(realization.num) || realization.num < 1) {
    stop("'realization.num' must be a positive integer")
  }

  if (length(daily.sim) != realization.num) {
    stop("Length of 'daily.sim' must equal 'realization.num'")
  }

  if (length(daily.obs) == 0) {
    stop("'daily.obs' must contain at least one grid cell")
  }

  if (!"precip" %in% variables) {
    stop("'variables' must include 'precip' for wet/dry spell analysis")
  }

  missing.vars <- setdiff(variables, names(daily.obs[[1]]))
  if (length(missing.vars) > 0) {
    stop("Variables not found in daily.obs: ", paste(missing.vars, collapse = ", "))
  }

  if (!"date" %in% names(daily.obs[[1]])) stop("'daily.obs' must contain a 'date' column")
  if (!"date" %in% names(daily.sim[[1]][[1]])) stop("'daily.sim' must contain a 'date' column")

  if (!is.numeric(wet.quantile) || wet.quantile <= 0 || wet.quantile >= 1) {
    stop("'wet.quantile' must be between 0 and 1")
  }

  if (!is.numeric(extreme.quantile) || extreme.quantile <= 0 || extreme.quantile >= 1) {
    stop("'extreme.quantile' must be between 0 and 1")
  }

  if (extreme.quantile <= wet.quantile) {
    stop("'extreme.quantile' must be greater than 'wet.quantile'")
  }

  invisible(TRUE)
}

#' Process observed weather data
#' @keywords internal
process_observed_data <- function(daily.obs, variables, n_grids,
                                  wet.quantile, extreme.quantile) {

  his.date <- daily.obs[[1]]$date

  his.datemat <- dplyr::tibble(
    date = his.date,
    year = as.integer(format(date, "%Y")),
    mon = as.integer(format(date, "%m")),
    day = as.integer(format(date, "%d")))

  his <- lapply(seq_len(n_grids), function(i) {
    df <- daily.obs[[i]][, variables, drop = FALSE]
    dplyr::bind_cols(his.datemat, df)
  }) %>%
    dplyr::bind_rows(.id = "id") %>%
    dplyr::mutate(id = as.integer(.data$id))

  mc.thresholds <- his %>%
    dplyr::group_by(.data$id, .data$mon) %>%
    dplyr::summarize(
      wet.th = {
        ppos <- .data$precip[is.finite(.data$precip) & .data$precip > 0]
        if (length(ppos) >= 5) stats::quantile(ppos, wet.quantile, names = FALSE, na.rm = TRUE) else 0
      },
      extreme.th = {
        ppos <- .data$precip[is.finite(.data$precip) & .data$precip > 0]
        if (length(ppos) >= 5) stats::quantile(ppos, extreme.quantile, names = FALSE, na.rm = TRUE) else 0
      },
      .groups = "drop")


  his.stats <- compute_timeseries_statistics(
    data = his,
    variables = variables,
    mc.thresholds = mc.thresholds
  )

  list(
    data = his,
    datemat = his.datemat,
    mc.thresholds = mc.thresholds,
    stats.season = his.stats$stats.season %>% dplyr::rename(Observed = .data$value),
    stats.mon.aavg = his.stats$stats.mon.aavg %>% dplyr::rename(Observed = .data$value),
    stats.annual.aavg = his.stats$stats.annual.aavg %>% dplyr::rename(Observed = .data$value),
    wetdry = his.stats$wetdry %>% dplyr::rename(Observed = .data$value),
    cor = his.stats$cor %>% dplyr::rename(Observed = .data$value),
    cor.cond = his.stats$cor.cond %>% dplyr::rename(Observed = .data$value)
  )
}

#' Process simulated weather data
#' @keywords internal
process_simulated_data <- function(daily.sim, realization.num, variables, mc.thresholds) {

  sim.datemat <- dplyr::tibble(
    date = daily.sim[[1]][[1]]$date,
    year = as.integer(format(date, "%Y")),
    mon = as.integer(format(date, "%m")),
    day = as.integer(format(date, "%d"))
  )

  sim <- lapply(seq_len(realization.num), function(i) {
    daily.sim[[i]] %>%
      dplyr::bind_rows(.id = "id") %>%
      dplyr::mutate(id = as.integer(.data$id)) %>%
      dplyr::left_join(sim.datemat, by = "date")
  })

  sim.stats.list <- lapply(seq_len(realization.num), function(i) {
    compute_timeseries_statistics(
      data = sim[[i]],
      variables = variables,
      mc.thresholds = mc.thresholds
    )
  })

  list(
    stats.season = dplyr::bind_rows(lapply(sim.stats.list, `[[`, "stats.season"), .id = "rlz") %>%
      dplyr::mutate(id = as.numeric(.data$id)) %>%
      dplyr::rename(Simulated = .data$value),

    stats.mon.aavg = dplyr::bind_rows(lapply(sim.stats.list, `[[`, "stats.mon.aavg"), .id = "rlz") %>%
      dplyr::rename(Simulated = .data$value),

    stats.annual.aavg = dplyr::bind_rows(lapply(sim.stats.list, `[[`, "stats.annual.aavg"), .id = "rlz") %>%
      dplyr::mutate(year = .data$year - min(.data$year) + 1) %>%
      dplyr::rename(Simulated = .data$value),

    cor = dplyr::bind_rows(lapply(sim.stats.list, `[[`, "cor"), .id = "rlz") %>%
      dplyr::rename(Simulated = .data$value),

    wetdry = dplyr::bind_rows(lapply(sim.stats.list, `[[`, "wetdry"), .id = "rlz") %>%
      dplyr::mutate(id = as.numeric(.data$id)) %>%
      dplyr::rename(Simulated = .data$value),

    cor.cond = dplyr::bind_rows(lapply(sim.stats.list, `[[`, "cor.cond"), .id = "rlz") %>%
      dplyr::rename(Simulated = .data$value)
  )
}

#' Compute fit metrics for each realization
#'
#' Computes Mean Absolute Error (MAE) for all key metrics.
#' MAE measures the average magnitude of errors between simulated and observed values.
#' Lower MAE indicates better fit.
#'
#' @param obs.results Observed data results from process_observed_data
#' @param sim.results Simulated data results from process_simulated_data
#' @param variables Character vector of variable names
#'
#' @keywords internal
compute_realization_fit_metrics <- function(obs.results, sim.results, variables) {

  metrics_list <- list()

  # --- 1. MAE for means (by variable, across all grid-months) ---
  mae_mean <- sim.results$stats.season %>%
    dplyr::filter(.data$stat == "mean") %>%
    dplyr::left_join(
      obs.results$stats.season %>% dplyr::select(.data$id, .data$mon, .data$variable, Observed),
      by = c("id", "mon", "variable")
    ) %>%
    dplyr::group_by(.data$rlz, .data$variable) %>%
    dplyr::summarize(
      mae_mean = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$mae_mean, names_prefix = "mae_mean_")

  # --- 2. MAE for SDs ---
  mae_sd <- sim.results$stats.season %>%
    dplyr::filter(.data$stat == "sd") %>%
    dplyr::left_join(
      obs.results$stats.season %>% dplyr::select(.data$id, .data$mon, .data$variable, Observed),
      by = c("id", "mon", "variable")
    ) %>%
    dplyr::group_by(.data$rlz, .data$variable) %>%
    dplyr::summarize(
      mae_sd = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$mae_sd, names_prefix = "mae_sd_")

  # --- 3. MAE for wet/dry day counts ---
  mae_wetdry_days <- sim.results$wetdry %>%
    dplyr::filter(.data$type == "days") %>%
    dplyr::left_join(
      obs.results$wetdry %>% dplyr::select(.data$id, .data$mon, .data$stat, Observed),
      by = c("id", "mon", "stat")
    ) %>%
    dplyr::group_by(.data$rlz, .data$stat) %>%
    dplyr::summarize(
      mae = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$stat, values_from = .data$mae, names_prefix = "mae_days_")

  # --- 4. MAE for spell lengths ---
  mae_spells <- sim.results$wetdry %>%
    dplyr::filter(.data$type == "spells") %>%
    dplyr::left_join(
      obs.results$wetdry %>% dplyr::select(.data$id, .data$mon, .data$stat, Observed),
      by = c("id", "mon", "stat")
    ) %>%
    dplyr::group_by(.data$rlz, .data$stat) %>%
    dplyr::summarize(
      mae = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$stat, values_from = .data$mae, names_prefix = "mae_spell_")

  # --- 5. MAE for cross-grid correlations ---
  mae_cor_crossgrid <- sim.results$cor %>%
    dplyr::filter(.data$variable1 == .data$variable2, .data$id1 != .data$id2) %>%
    dplyr::left_join(
      obs.results$cor %>%
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
  mae_cor_intervariable <- sim.results$cor %>%
    dplyr::filter(.data$id1 == .data$id2, .data$variable1 != .data$variable2) %>%
    dplyr::left_join(
      obs.results$cor %>%
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


#' Compute time series statistics
#' @keywords internal
compute_timeseries_statistics <- function(data, variables, mc.thresholds) {

  stat.funs <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)
  year.num <- length(unique(format(data$date, "%Y")))

  stats.season <- compute_grouped_statistics(
    df = data,
    variables = variables,
    group.vars = c("id", "mon"),
    stat.funs = stat.funs
  )

  stats.mon.aavg <- compute_grouped_statistics(
    df = data,
    variables = variables,
    group.vars = c("year", "mon"),
    stat.funs = stat.funs
  )

  stats.annual.aavg <- compute_grouped_statistics(
    df = data,
    variables = variables,
    group.vars = "year",
    stat.funs = stat.funs
  ) %>%
    dplyr::mutate(year = .data$year - min(.data$year) + 1)

  wetdry <- data %>%
    dplyr::left_join(mc.thresholds, by = c("id", "mon")) %>%
    dplyr::group_by(.data$id, .data$mon) %>%
    dplyr::summarize(
      Wet_days = {
        th <- .data$wet.th[1]
        # If threshold is 0 (fallback or quantile result), treat wet as strictly > 0
        if (!is.finite(th) || th <= 0) sum(.data$precip > 0, na.rm = TRUE) / year.num
        else sum(.data$precip > th, na.rm = TRUE) / year.num
      },
      Dry_days = {
        th <- .data$wet.th[1]
        if (!is.finite(th) || th <= 0) sum(.data$precip <= 0, na.rm = TRUE) / year.num
        else sum(.data$precip <= th, na.rm = TRUE) / year.num
      },
      Dry_spells = {
        th <- .data$wet.th[1]
        if (!is.finite(th) || th <= 0) mean_spell_length(.data$precip, threshold = 0, below = TRUE)
        else mean_spell_length(.data$precip, threshold = th, below = TRUE)
      },
      Wet_spells = {
        th <- .data$wet.th[1]
        if (!is.finite(th) || th <= 0) mean_spell_length(.data$precip, threshold = 0, below = FALSE)
        else mean_spell_length(.data$precip, threshold = th, below = FALSE)
      },
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(.data$Wet_days, .data$Dry_days, .data$Dry_spells, .data$Wet_spells),
      names_to = "stat.full",
      values_to = "value"
    ) %>%
    tidyr::separate(.data$stat.full, into = c("stat", "type"), sep = "_") %>%
    dplyr::mutate(variable = "precip", .after = .data$mon)


  cor.data <- compute_correlation_matrix_anom(data, variables)

  cor.cond <- compute_conditional_precip_cor(
    data = data,
    variables = variables,
    mc.thresholds = mc.thresholds,
    wet_def = "monthly_quantile",
    use_anom = TRUE,
    use_log_precip_on_wet = TRUE,
    method = "pearson",
    min_pairs = 50
  )

  list(
    stats.season = stats.season,
    stats.mon.aavg = stats.mon.aavg,
    stats.annual.aavg = stats.annual.aavg,
    wetdry = wetdry,
    cor = cor.data,
    cor.cond = cor.cond
  )
}

#' Compute grouped statistics
#' @keywords internal
compute_grouped_statistics <- function(df, variables, group.vars, stat.funs) {

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group.vars))) %>%
    dplyr::summarize(
      dplyr::across(dplyr::all_of(variables), stat.funs, .names = "{.col}:{.fn}"),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(group.vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    tidyr::separate(.data$variable, into = c("variable", "stat"), sep = ":") %>%
    dplyr::mutate(
      stat = factor(.data$stat, levels = names(stat.funs)),
      variable = as.character(.data$variable)
    )
}

#' Compute correlation matrix on monthly anomalies
#' @keywords internal
compute_correlation_matrix_anom <- function(data, variables) {

  dat2 <- data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(variables),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::group_by(.data$id, .data$mon, .data$variable) %>%
    dplyr::mutate(value = value - mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    tidyr::unite("id.variable", .data$id, .data$variable, sep = ":") %>%
    dplyr::select(.data$date, .data$id.variable, .data$value) %>%
    tidyr::pivot_wider(names_from = .data$id.variable, values_from = .data$value) %>%
    dplyr::arrange(.data$date)

  mat <- as.matrix(dat2 %>% dplyr::select(-.data$date))
  cmat <- stats::cor(mat, use = "pairwise.complete.obs", method = "pearson")

  tri <- upper.tri(cmat, diag = FALSE)

  dplyr::tibble(
    id.variable1 = colnames(cmat)[row(cmat)[tri]],
    id.variable2 = colnames(cmat)[col(cmat)[tri]],
    value = cmat[tri]
  ) %>%
    tidyr::separate(.data$id.variable1, c("id1", "variable1"), sep = ":") %>%
    tidyr::separate(.data$id.variable2, c("id2", "variable2"), sep = ":")
}

#' Conditional precip-X correlations within each grid
#' @keywords internal
compute_conditional_precip_cor <- function(data,
                                           variables,
                                           mc.thresholds = NULL,
                                           wet_def = c("gt0", "monthly_quantile"),
                                           use_anom = TRUE,
                                           use_log_precip_on_wet = TRUE,
                                           method = "pearson",
                                           min_pairs = 50) {

  wet_def <- match.arg(wet_def)

  # keep only variables we actually have
  if (!("precip" %in% variables)) stop("variables must include 'precip'", call. = FALSE)
  targets <- setdiff(variables, "precip")
  if (length(targets) == 0) {
    return(dplyr::tibble(
      id1 = integer(), variable1 = character(),
      id2 = integer(), variable2 = character(),
      regime = character(), transform = character(),
      value = numeric(), n = integer()
    ))
  }

  # Optionally de-seasonalize (monthly mean anomalies) for all variables
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

  # Wide format back (per grid) to compute pairwise cor easily with filtering
  dat_wide <- dat %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$value)

  # Join wet thresholds if needed
  if (wet_def == "monthly_quantile") {
    if (is.null(mc.thresholds)) stop("mc.thresholds required when wet_def='monthly_quantile'", call. = FALSE)
    dat_wide <- dat_wide %>% dplyr::left_join(mc.thresholds, by = c("id", "mon"))
    wet_flag <- dat_wide$precip >= dat_wide$wet.th
  } else {
    wet_flag <- dat_wide$precip > 0
  }
  dry_flag <- !wet_flag & !is.na(dat_wide$precip)

  # Build results
  out <- lapply(targets, function(v) {

    # all days (raw/anom precip)
    res_all <- .safe_cor(dat_wide$precip, dat_wide[[v]])

    # wet days (optionally log1p on precip)
    if (use_log_precip_on_wet) {
      p_wet <- log1p(pmax(dat_wide$precip, 0))
      res_wet <- .safe_cor(p_wet[wet_flag], dat_wide[[v]][wet_flag])
      tr_wet <- "log1p_precip"
    } else {
      res_wet <- .safe_cor(dat_wide$precip[wet_flag], dat_wide[[v]][wet_flag])
      tr_wet <- "precip"
    }

    # dry days (precip is ~0; correlation often meaningless for precip itself, but can be informative if precip anomalies exist)
    res_dry <- .safe_cor(dat_wide$precip[dry_flag], dat_wide[[v]][dry_flag])

    dplyr::tibble(
      id1 = dat_wide$id[1],  # placeholder, fixed below
      variable1 = "precip",
      id2 = dat_wide$id[1],  # placeholder, fixed below
      variable2 = v,
      regime = c("all", "wet", "dry"),
      transform = c("precip", tr_wet, "precip"),
      value = c(res_all["value"], res_wet["value"], res_dry["value"]),
      n = c(as.integer(res_all["n"]), as.integer(res_wet["n"]), as.integer(res_dry["n"]))
    )
  }) %>%
    dplyr::bind_rows()

  # Fix id columns per-grid: compute per id properly
  # We need per-grid computation, not collapsing all ids.
  # So do it by splitting dat_wide by id:
  by_id <- split(dat_wide, dat_wide$id)

  out2 <- lapply(by_id, function(df_id) {

    if (wet_def == "monthly_quantile") {
      wet_flag <- df_id$precip >= df_id$wet.th
    } else {
      wet_flag <- df_id$precip > 0
    }
    dry_flag <- !wet_flag & !is.na(df_id$precip)

    .safe_cor <- function(x, y) {
      ok <- is.finite(x) & is.finite(y)
      n <- sum(ok)
      if (n < min_pairs) return(c(value = NA_real_, n = n))
      c(value = stats::cor(x[ok], y[ok], method = method), n = n)
    }

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

  out2
}


#' Generate symmetric dummy points for equal axis scaling
#' @keywords internal
generate_symmetric_dummy_points <- function(df, facet.var = "variable",
                                            x.col = "Observed", y.col = "Simulated") {

  df %>%
    dplyr::filter(
      !is.na(.data[[x.col]]),
      !is.na(.data[[y.col]]),
      is.finite(.data[[x.col]]),
      is.finite(.data[[y.col]])
    ) %>%
    dplyr::mutate(
      max.val = pmax(.data[[x.col]], .data[[y.col]]),
      min.val = pmin(.data[[x.col]], .data[[y.col]])
    ) %>%
    dplyr::group_by(.data[[facet.var]]) %>%
    dplyr::summarize(
      minlim = min(.data$min.val, na.rm = TRUE),
      maxlim = max(.data$max.val, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::expand_grid(lim = c("min", "max")) %>%
    dplyr::mutate(
      Observed = ifelse(.data$lim == "min", .data$minlim, .data$maxlim),
      Simulated = ifelse(.data$lim == "min", .data$minlim, .data$maxlim)
    ) %>%
    dplyr::select(dplyr::all_of(facet.var), .data$Observed, .data$Simulated)
}


#' Filter correlations by type (cross-grid or inter-variable)
#' @keywords internal
filter_correlations <- function(cor_data, same_var, same_id, allowed_pairs, pair_type) {

  result <- cor_data %>%
    dplyr::filter(
      (.data$variable1 == .data$variable2) == same_var,
      (.data$id1 == .data$id2) == same_id
    )

  if (pair_type == "id") {
    result <- result %>%
      dplyr::mutate(
        pair = paste(pmin(.data$id1, .data$id2), pmax(.data$id1, .data$id2), sep = ":")
      ) %>%
      dplyr::filter(.data$pair %in% allowed_pairs) %>%
      dplyr::select(-.data$pair)

  } else if (pair_type == "variable") {
    result <- result %>%
      dplyr::mutate(
        pair = paste(pmin(.data$variable1, .data$variable2),
                     pmax(.data$variable1, .data$variable2), sep = ":"),
        variable = .data$pair
      ) %>%
      dplyr::filter(.data$pair %in% allowed_pairs) %>%
      dplyr::select(-.data$pair)
  }

  result
}

#' Prepare data for plotting
#' @keywords internal
prepare_plot_data <- function(obs.results, sim.results, variables, n_grids) {

  # Define statistical functions for reference
  stat.funs <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)

  # Merge seasonal statistics
  daily.stats.season <- sim.results$stats.season %>%
    dplyr::left_join(
      obs.results$stats.season,
      by = c("id", "mon", "variable", "stat")
    ) %>%
    dplyr::mutate(
      variable = factor(.data$variable),
      stat = factor(.data$stat, levels = names(stat.funs))
    )

  # Merge correlations (unconditional anomalies)
  stats.allcor <- sim.results$cor %>%
    dplyr::left_join(
      obs.results$cor,
      by = c("id1", "variable1", "id2", "variable2")
    )

  # Merge conditional precip correlations (within-grid; new)
  stats.allcor.cond <- sim.results$cor.cond %>%
    dplyr::left_join(
      obs.results$cor.cond,
      by = c("id1", "variable1", "id2", "variable2", "regime", "transform")
    )

  # ---------------------------------------------------------------------------
  # Robust ID + variable pairing without assuming ids are 1..n
  # ---------------------------------------------------------------------------

  # Prefer IDs from the observed long data if available; otherwise infer from correlations
  ids_from_obs <- NULL
  if (!is.null(obs.results$data) && "id" %in% names(obs.results$data)) {
    ids_from_obs <- sort(unique(obs.results$data$id))
  }

  ids_from_cor <- NULL
  if (nrow(stats.allcor) > 0) {
    ids_from_cor <- sort(unique(c(stats.allcor$id1, stats.allcor$id2)))
  }

  ids <- ids_from_obs
  if (is.null(ids) || length(ids) == 0) ids <- ids_from_cor
  if (is.null(ids) || length(ids) < 2) {
    # no meaningful cross-grid pairs can be formed
    ids <- unique(ids)
  }

  # Allowed grid-pair keys (order-invariant)
  id_pairs_allowed <- character(0)
  if (!is.null(ids) && length(ids) >= 2) {
    id_pairs_allowed <- apply(utils::combn(ids, 2), 2, function(x) paste(x[1], x[2], sep = ":"))
  }

  # Allowed variable-pair keys (order-invariant)
  vars_use <- unique(as.character(variables))
  var_pairs_allowed <- character(0)
  if (length(vars_use) >= 2) {
    var_pairs_allowed <- apply(utils::combn(sort(vars_use), 2), 2, function(x) paste(x[1], x[2], sep = ":"))
  }

  # ---------------------------------------------------------------------------
  # Cross-grid correlations: same variable, different grids
  # ---------------------------------------------------------------------------

  stats.crosscor <- filter_correlations(cor_data = stats.allcor,
                                        same_var = TRUE, same_id = FALSE, allowed_pairs = id_pairs_allowed,
                                        pair_type = "id")

  # ---------------------------------------------------------------------------
  # Inter-variable correlations: same grid, different variables
  # ---------------------------------------------------------------------------

  stats.intercor <- filter_correlations(cor_data = stats.allcor,
                                        same_var = FALSE, same_id = TRUE, allowed_pairs = var_pairs_allowed,
                                        pair_type = "variable")

  # ---------------------------------------------------------------------------
  # Wet/dry statistics
  # ---------------------------------------------------------------------------

  stats.wetdry <- sim.results$wetdry %>%
    dplyr::left_join(
      obs.results$wetdry,
      by = c("id", "mon", "variable", "stat", "type")
    ) %>%
    dplyr::mutate(stat = factor(.data$stat, levels = c("Dry", "Wet")))

  list(
    daily.stats.season = daily.stats.season,
    stats.mon.aavg.sim = sim.results$stats.mon.aavg,
    stats.mon.aavg.obs = obs.results$stats.mon.aavg,
    stats.annual.aavg.sim = sim.results$stats.annual.aavg,
    stats.annual.aavg.obs = obs.results$stats.annual.aavg,
    stats.crosscor = stats.crosscor,
    stats.intercor = stats.intercor,
    stats.wetdry = stats.wetdry,
    stats.precip_cor_cond = stats.allcor.cond
  )
}
