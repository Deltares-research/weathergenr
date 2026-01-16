#' Evaluate Stochastic Weather Generator Performance
#'
#' @description
#' Run a comprehensive diagnostic evaluation comparing synthetic weather
#' simulations against historical observations across multiple grid cells.
#' Computes summary statistics, wet/dry day counts, spell lengths, and
#' inter-site correlations, and generates diagnostic plots to assess
#' stochastic weather generator performance.
#'
#' @details
#' The function standardizes simulated and observed series to full-year windows
#' (after leap-day removal), optionally subsamples grid cells to manage memory,
#' and returns diagnostic plots alongside a summarized fit table. Missing values
#' are ignored in summary statistics (`na.rm = TRUE`), and correlations use
#' pairwise complete observations. Use `seed` for reproducible subsampling and
#' window selection; the original RNG state is restored on exit.
#'
#' @param daily_sim List of simulated weather realizations. Each element should be
#'   a list of data frames (one per grid cell), containing daily values and a `date` column.
#' @param daily_obs List of observed weather data frames (one per grid cell). Each
#'   should contain a `date` column and the variables specified in `variables`.
#' @param variables Character vector of variable names to evaluate (e.g., `c("precip", "temp")`).
#'   Must include "precip" for wet/dry spell analysis.
#' @param variable_labels Optional character vector of variable labels for plots.
#'   Defaults to `variables` if `NULL`.
#' @param n_realizations Integer. Number of synthetic realizations in `daily_sim`.
#' @param wet_quantile Numeric between 0 and 1. Quantile threshold for wet days
#'   (default = 0.2).
#' @param extreme_quantile Numeric between 0 and 1. Quantile threshold for extremely
#'   wet days (default = 0.8).
#' @param output_path Character. Directory path to save generated plots. If `NULL`,
#'   plots are not saved to disk.
#' @param save_plots Logical. Whether to save plots to `output_path` (default = `TRUE`).
#' @param show_title Logical. Whether to display titles in plots (default = `TRUE`).
#' @param verbose Logical. Whether to emit console messages and the fit summary table.
#' @param max_grids Integer. Maximum number of grid cells to evaluate (default = 25).
#'   If more grids are provided, a random subsample is used to control memory.
#' @param seed Optional integer. Random seed for reproducible grid subsampling and
#'   year window selection. If NULL, results will vary between runs.
#'
#' @return A named list of `ggplot2` plot objects with class "weather_assessment".
#'   The returned object also contains attributes:
#'   \itemize{
#'     \item \code{fit_summary}: data frame of per-realization fit metrics and ranks.
#'     \item \code{metadata}: list with \code{n_grids}, \code{n_realizations},
#'       \code{variables}, and \code{assessment_date}.
#'   }
#'
#' @examples
#' set.seed(1)
#' dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
#' dates <- dates[format(dates, "%m-%d") != "02-29"]
#' obs_grid <- list(data.frame(
#'   date = dates,
#'   precip = rgamma(length(dates), shape = 2, scale = 2),
#'   temp = rnorm(length(dates), mean = 10, sd = 3)
#' ))
#' sim_grid <- list(list(data.frame(
#'   date = dates,
#'   precip = rgamma(length(dates), shape = 2, scale = 2),
#'   temp = rnorm(length(dates), mean = 10, sd = 3)
#' )))
#' out <- evaluate_weather_generator(
#'   daily_sim = sim_grid,
#'   daily_obs = obs_grid,
#'   variables = c("precip", "temp"),
#'   n_realizations = 1,
#'   output_path = NULL,
#'   save_plots = FALSE,
#'   show_title = FALSE
#' )
#' class(out)
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom e1071 skewness
#' @importFrom stats cor sd quantile
#' @importFrom utils combn
#' @export
evaluate_weather_generator <- function(
    daily_sim = NULL,
    daily_obs = NULL,
    variables = NULL,
    variable_labels = NULL,
    n_realizations = NULL,
    wet_quantile = 0.2,
    extreme_quantile = 0.8,
    output_path = NULL,
    save_plots = TRUE,
    show_title = TRUE,
    verbose = TRUE,
    max_grids = 25,
    seed = NULL
) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  .validate_assessment_inputs(
    daily_sim = daily_sim,
    daily_obs = daily_obs,
    variables = variables,
    n_realizations = n_realizations,
    wet_quantile = wet_quantile,
    extreme_quantile = extreme_quantile,
    verbose = verbose
  )

  if (!.is_int_scalar(max_grids) || max_grids < 1L) {
    stop("'max_grids' must be a positive integer", call. = FALSE)
  }

  # ============================================================================
  # RNG STATE MANAGEMENT
  # ============================================================================

  if (!is.null(seed)) {
    if (!.is_int_scalar(seed)) {
      stop("'seed' must be NULL or a single integer", call. = FALSE)
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

  .log(
    paste0(
      "[VALIDATE] Start | grids = {format(length(daily_obs), big.mark = ',')} | ",
      "realizations = {format(n_realizations, big.mark = ',')}"
    ),
    verbose = verbose
  )
  .log(
    paste0(
      "[VALIDATE] Variables = {paste(variables, collapse = ',')}"
    ),
    verbose = verbose
  )
  .log(
    paste0("[VALIDATE] Parameters: wet.q = {wet_quantile} | extreme.q = {extreme_quantile}"),
    verbose = verbose
  )


  if (is.null(variable_labels)) variable_labels <- variables

  if (!is.null(output_path)) {
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    save_plots <- FALSE
  }

  options(dplyr.summarise.inform = FALSE, tidyverse.quiet = TRUE)

  plot_config <- list(
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

  grid_count <- length(daily_obs)
  grid_count_original <- grid_count

  if (grid_count > max_grids) {

    sel.grids <- sort(sample(seq_len(grid_count), max_grids))

    daily_obs <- daily_obs[sel.grids]
    daily_sim <- lapply(daily_sim, function(rlz) rlz[sel.grids])

    grid_count <- length(daily_obs)

    if (isTRUE(verbose) && requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn(
        "[VALIDATE] Grid count reduced from {format(grid_count_original, big.mark = ',')} to {format(grid_count, big.mark = ',')} for memory control."
      )
    }
  }


  # ============================================================================
  # STANDARDIZE PERIODS (FULL YEARS + MATCH LENGTH VIA RANDOM WINDOW)
  # ============================================================================

  .log(
    "[VALIDATE] Standardizing obs/sim periods to full years and equal length",
    verbose = verbose
  )

  std <- .align_obs_sim_periods(
    daily_obs = daily_obs,
    daily_sim = daily_sim,
    n_realizations = n_realizations,
    variables = variables
  )

  daily_obs <- std$daily_obs
  daily_sim <- std$daily_sim

  .log(
    paste0(
      "[VALIDATE] Standardized period | ",
      "Obs = ", std$obs_year_start, "-", std$obs_year_end, " | ",
      "Sim = ", std$sim_year_start, "-", std$sim_year_end
    ),
    verbose = verbose
  )

  # ============================================================================
  # PROCESS OBSERVED DATA
  # ============================================================================

  .log("[VALIDATE] Processing observed data", verbose = verbose)

  obs_results <- .summarize_observed_data(
    daily_obs = daily_obs,
    variables = variables,
    grid_count = grid_count,
    wet_quantile = wet_quantile,
    extreme_quantile = extreme_quantile
  )

  # ============================================================================
  # PROCESS SIMULATED DATA
  # ============================================================================

  .log(
    "[VALIDATE] Processing simulated data ({format(n_realizations, big.mark = ',')} realizations)",
    verbose = verbose
  )

  sim_results <- .summarize_simulated_data(
    daily_sim = daily_sim,
    n_realizations = n_realizations,
    variables = variables,
    mc_thresholds = obs_results$mc_thresholds
  )

  # ============================================================================
  # MERGE AND PREPARE PLOT DATA
  # ============================================================================

  .log("[VALIDATE] Preparing diagnostic data for plotting", verbose = verbose)

  plot_data <- .build_plot_data(
    obs_results = obs_results,
    sim_results = sim_results,
    variables = variables
  )

  # ============================================================================
  # GENERATE DIAGNOSTIC PLOTS
  # ============================================================================

  .log("[VALIDATE] Generating diagnostic plots", verbose = verbose)

  plots <- create_all_diagnostic_plots(
    plot_data = plot_data,
    plot_config = plot_config,
    variables = variables,
    show_title = show_title,
    save_plots = save_plots,
    output_path = output_path
  )

  .log(
    "[VALIDATE] Generated {format(length(plots), big.mark = ',')} diagnostic plots.",
    verbose = verbose
  )
  if (save_plots) {
    .log("[VALIDATE] Plots saved to: {output_path}", verbose = verbose)
  }

  # ============================================================================
  # COMPUTE FIT METRICS SUMMARY TABLE
  # ============================================================================

  .log("[VALIDATE] Computing fit metrics for all realizations", verbose = verbose)

  fit_summary <- .summarize_realization_fit(
    obs_results = obs_results,
    sim_results = sim_results,
    variables = variables
  )

  # ============================================================================
  # DISPLAY FIT SUMMARY TABLE
  # ============================================================================

  .log("[VALIDATE] Displaying fit assessment summary", verbose = verbose)

  if (isTRUE(verbose)) {
    .print_fit_summary_table(fit_summary)
  }

  .log("[VALIDATE] Assessment completed successfully", verbose = verbose)

  structure(
    plots,
    class = c("weather_assessment", "list"),
    fit_summary = fit_summary,
    metadata = list(
      n_grids = grid_count,
      n_realizations = n_realizations,
      variables = variables,
      assessment_date = Sys.Date()
    )
  )
}

# ==============================================================================
# FIT SUMMARY DISPLAY
# ==============================================================================

#' Display fit summary table for all realizations
#'
#' Formats and prints a compact summary of fit metrics across realizations.
#' This is primarily intended for interactive review of model performance.
#'
#' @param fit_summary Data frame returned by `.summarize_realization_fit()`.
#' @return Invisibly returns `fit_summary`.
#' @keywords internal
.print_fit_summary_table <- function(fit_summary) {

  if (is.null(fit_summary) || nrow(fit_summary) == 0) {
    cat("\n")
    cat("Warning: no fit metrics computed.\n")
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
    vapply(
      col_names,
      function(nm) {
        if (nm %in% c("rlz", "rank")) {
          max(nchar(format(display_data[[nm]], big.mark = ",")))
        } else {
          max(nchar(format(display_data[[nm]], digits = 4, nsmall = 4)))
        }
      },
      integer(1)
    )
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
        row_values[j] <- format(
          as.integer(display_data[i, j]),
          width = col_widths[j],
          justify = "right",
          big.mark = ","
        )
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
  cat("  - Best realization  : ", format(best_rlz, big.mark = ","), " (score = ",
      format(round(fit_summary$overall_score[1], 4), nsmall = 4), ")\n", sep = "")
  cat("  - Worst realization : ", format(worst_rlz, big.mark = ","), " (score = ",
      format(round(fit_summary$overall_score[nrow(fit_summary)], 4), nsmall = 4), ")\n", sep = "")
  cat("  - Median score      : ", format(round(median_score, 4), nsmall = 4), "\n", sep = "")
  cat("\n")

  invisible(fit_summary)
}

# ==============================================================================
# HELPERS
# ==============================================================================

#' Find full (365-day) years after leap-day removal and return longest contiguous block
#'
#' @param dates Date vector (leap days should already be removed).
#' @return List with a `years` integer vector. Empty if no full-year blocks exist.
#' @keywords internal
.get_full_year_run <- function(dates) {

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

#' Pick a random contiguous window of years of length `window_years` from a year vector
#'
#' @param years Integer vector of available years.
#' @param window_years Integer number of contiguous years to sample.
#' @return Integer vector of selected years in chronological order.
#' @keywords internal
.pick_contiguous_year_window <- function(years, window_years) {

  years <- as.integer(years)
  if (length(years) < window_years) stop("Not enough full years to pick window.", call. = FALSE)
  if (window_years == length(years)) return(years)

  start_idx <- sample.int(length(years) - window_years + 1L, size = 1L)
  years[start_idx:(start_idx + window_years - 1L)]
}

#' Filter a single grid data frame to selected years (keeps original columns)
#'
#' @param df Data frame with a `date` column.
#' @param years_keep Integer vector of years to retain.
#' @return Filtered data frame with the same columns as `df`.
#' @keywords internal
.filter_grid_years <- function(df, years_keep) {

  if (!("date" %in% names(df))) stop("df must contain 'date'", call. = FALSE)

  yrs <- as.integer(format(df$date, "%Y"))
  df[yrs %in% years_keep, , drop = FALSE]
}

#' Standardize obs and sim to the same full-year length via random windowing
#'
#' Applies the chosen window consistently across all grids and all realizations
#' after removing leap days. The window is drawn from the longest contiguous
#' full-year blocks available in each series.
#'
#' @param daily_obs List of observed data frames (one per grid).
#' @param daily_sim List of simulated realizations; each realization is a list of grids.
#' @param n_realizations Integer number of realizations.
#' @param variables Character vector of variables used for evaluation.
#' @return List with standardized `daily_obs`, `daily_sim`, and window metadata.
#' @keywords internal
.align_obs_sim_periods <- function(daily_obs, daily_sim, n_realizations, variables) {

  # 1) Remove leap days consistently (based on each series' date vector)
  #    Assumption: all grids share the same date vector within obs and within each realization.
  obs_dates <- daily_obs[[1]]$date
  sim_dates <- daily_sim[[1]][[1]]$date

  obs_leap <- find_leap_day_indices(obs_dates)
  sim_leap <- find_leap_day_indices(sim_dates)

  if (!is.null(obs_leap)) {
    daily_obs <- lapply(daily_obs, function(df) df[-obs_leap, , drop = FALSE])
    obs_dates <- daily_obs[[1]]$date
  }

  if (!is.null(sim_leap)) {
    daily_sim <- lapply(daily_sim, function(rlz) {
      lapply(rlz, function(df) df[-sim_leap, , drop = FALSE])
    })
    sim_dates <- daily_sim[[1]][[1]]$date
  }

  # 2) Identify longest contiguous full-year blocks
  obs_block <- .get_full_year_run(obs_dates)$years
  sim_block <- .get_full_year_run(sim_dates)$years

  if (length(obs_block) == 0) stop("Observed series has no complete 365-day years after leap-day removal.", call. = FALSE)
  if (length(sim_block) == 0) stop("Simulated series has no complete 365-day years after leap-day removal.", call. = FALSE)

  n_years <- min(length(obs_block), length(sim_block))

  # 3) Pick windows (randomly within each block if longer than needed)
  obs_years_keep <- .pick_contiguous_year_window(obs_block, window_years = n_years)
  sim_years_keep <- .pick_contiguous_year_window(sim_block, window_years = n_years)

  # 4) Apply year filtering
  daily_obs2 <- lapply(daily_obs, .filter_grid_years, years_keep = obs_years_keep)

  daily_sim2 <- lapply(seq_len(n_realizations), function(i) {
    lapply(daily_sim[[i]], .filter_grid_years, years_keep = sim_years_keep)
  })

  # Basic integrity check: after filtering, all grids should have same nrow within obs and sim
  obs_n <- nrow(daily_obs2[[1]])
  sim_n <- nrow(daily_sim2[[1]][[1]])

  if (obs_n != sim_n) {
    stop(
      "Period standardization failed: obs and sim window lengths differ after filtering. ",
      "obs_n=", obs_n, " sim_n=", sim_n,
      call. = FALSE
    )
  }

  list(
    daily_obs = daily_obs2,
    daily_sim = daily_sim2,
    n_years = n_years,
    obs_year_start = min(obs_years_keep),
    obs_year_end = max(obs_years_keep),
    sim_year_start = min(sim_years_keep),
    sim_year_end = max(sim_years_keep)
  )
}

#' Validate inputs for weather assessment
#' @keywords internal
.validate_assessment_inputs <- function(
    daily_sim,
    daily_obs,
    variables,
    n_realizations,
    wet_quantile,
    extreme_quantile,
    verbose
) {

  if (is.null(daily_sim)) stop("'daily_sim' must not be NULL", call. = FALSE)
  if (is.null(daily_obs)) stop("'daily_obs' must not be NULL", call. = FALSE)
  if (is.null(variables)) stop("'variables' must not be NULL", call. = FALSE)
  if (is.null(n_realizations)) stop("'n_realizations' must not be NULL", call. = FALSE)

  if (!is.list(daily_sim)) stop("'daily_sim' must be a list", call. = FALSE)
  if (!is.list(daily_obs)) stop("'daily_obs' must be a list", call. = FALSE)
  if (!is.character(variables)) stop("'variables' must be a character vector", call. = FALSE)
  if (!.is_int_scalar(n_realizations) || n_realizations < 1L) {
    stop("'n_realizations' must be a positive integer", call. = FALSE)
  }

  if (length(daily_sim) != n_realizations) {
    stop("Length of 'daily_sim' must equal 'n_realizations'", call. = FALSE)
  }

  if (length(daily_obs) == 0) {
    stop("'daily_obs' must contain at least one grid cell", call. = FALSE)
  }

  if (!"precip" %in% variables) {
    stop("'variables' must include 'precip' for wet/dry spell analysis", call. = FALSE)
  }

  missing.vars <- setdiff(variables, names(daily_obs[[1]]))
  if (length(missing.vars) > 0) {
    stop("Variables not found in daily_obs: ", paste(missing.vars, collapse = ", "), call. = FALSE)
  }

  if (!"date" %in% names(daily_obs[[1]])) {
    stop("'daily_obs' must contain a 'date' column", call. = FALSE)
  }
  if (!"date" %in% names(daily_sim[[1]][[1]])) {
    stop("'daily_sim' must contain a 'date' column", call. = FALSE)
  }

  if (!is.numeric(wet_quantile) || !is.finite(wet_quantile) ||
      wet_quantile <= 0 || wet_quantile >= 1) {
    stop("'wet_quantile' must be between 0 and 1", call. = FALSE)
  }

  if (!is.numeric(extreme_quantile) || !is.finite(extreme_quantile) ||
      extreme_quantile <= 0 || extreme_quantile >= 1) {
    stop("'extreme_quantile' must be between 0 and 1", call. = FALSE)
  }

  if (extreme_quantile <= wet_quantile) {
    stop("'extreme_quantile' must be greater than 'wet_quantile'", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE or FALSE", call. = FALSE)
  }

  invisible(TRUE)
}

#' Summarize observed weather data
#' @keywords internal
.summarize_observed_data <- function(daily_obs, variables, grid_count,
                                    wet_quantile, extreme_quantile) {

  his.date <- daily_obs[[1]]$date

  his.datemat <- dplyr::tibble(
    date = his.date,
    year = as.integer(format(date, "%Y")),
    mon = as.integer(format(date, "%m")),
    day = as.integer(format(date, "%d")))

  his <- lapply(seq_len(grid_count), function(i) {
    df <- daily_obs[[i]][, variables, drop = FALSE]
    dplyr::bind_cols(his.datemat, df)
  }) %>%
    dplyr::bind_rows(.id = "id") %>%
    dplyr::mutate(id = as.integer(.data$id))

  mc_thresholds <- his %>%
    dplyr::group_by(.data$id, .data$mon) %>%
    dplyr::summarize(
      wet.th = {
        ppos <- .data$precip[is.finite(.data$precip) & .data$precip > 0]
        if (length(ppos) >= 5) stats::quantile(ppos, wet_quantile, names = FALSE, na.rm = TRUE) else 0
      },
      extreme.th = {
        ppos <- .data$precip[is.finite(.data$precip) & .data$precip > 0]
        if (length(ppos) >= 5) stats::quantile(ppos, extreme_quantile, names = FALSE, na.rm = TRUE) else 0
      },
      .groups = "drop")


  his.stats <- .compute_timeseries_stats(
    data = his,
    variables = variables,
    mc_thresholds = mc_thresholds
  )

  list(
    data = his,
    datemat = his.datemat,
    mc_thresholds = mc_thresholds,
    stats.season = his.stats$stats.season %>% dplyr::rename(Observed = .data$value),
    stats.mon.aavg = his.stats$stats.mon.aavg %>% dplyr::rename(Observed = .data$value),
    stats.annual.aavg = his.stats$stats.annual.aavg %>% dplyr::rename(Observed = .data$value),
    wetdry = his.stats$wetdry %>% dplyr::rename(Observed = .data$value),
    cor = his.stats$cor %>% dplyr::rename(Observed = .data$value),
    cor.cond = his.stats$cor.cond %>% dplyr::rename(Observed = .data$value)
  )
}

#' Summarize simulated weather data
#' @keywords internal
.summarize_simulated_data <- function(daily_sim, n_realizations, variables, mc_thresholds) {

  sim.datemat <- dplyr::tibble(
    date = daily_sim[[1]][[1]]$date,
    year = as.integer(format(date, "%Y")),
    mon = as.integer(format(date, "%m")),
    day = as.integer(format(date, "%d"))
  )

  sim <- lapply(seq_len(n_realizations), function(i) {
    daily_sim[[i]] %>%
      dplyr::bind_rows(.id = "id") %>%
      dplyr::mutate(id = as.integer(.data$id)) %>%
      dplyr::left_join(sim.datemat, by = "date")
  })

  sim.stats.list <- lapply(seq_len(n_realizations), function(i) {
    .compute_timeseries_stats(
      data = sim[[i]],
      variables = variables,
      mc_thresholds = mc_thresholds
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
#' @param obs_results Observed data results from .summarize_observed_data
#' @param sim_results Simulated data results from .summarize_simulated_data
#' @param variables Character vector of variable names
#'
#' @keywords internal
.summarize_realization_fit <- function(obs_results, sim_results, variables) {

  metrics_list <- list()

  # --- 1. MAE for means (by variable, across all grid-months) ---
  mae_mean <- sim_results$stats.season %>%
    dplyr::filter(.data$stat == "mean") %>%
    dplyr::left_join(
      obs_results$stats.season %>% dplyr::select(.data$id, .data$mon, .data$variable, Observed),
      by = c("id", "mon", "variable")
    ) %>%
    dplyr::group_by(.data$rlz, .data$variable) %>%
    dplyr::summarize(
      mae_mean = mean(abs(.data$Simulated - .data$Observed), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$mae_mean, names_prefix = "mae_mean_")

  # --- 2. MAE for SDs ---
  mae_sd <- sim_results$stats.season %>%
    dplyr::filter(.data$stat == "sd") %>%
    dplyr::left_join(
      obs_results$stats.season %>% dplyr::select(.data$id, .data$mon, .data$variable, Observed),
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


#' Compute time series statistics
#'
#' @param data Data frame with daily values and a `date` column.
#' @param variables Character vector of variable names to summarize.
#' @param mc_thresholds Data frame of wet-day thresholds by grid and month.
#' @return List of seasonal, monthly, and annual stats plus wet/dry and correlation data.
#' @keywords internal
.compute_timeseries_stats <- function(data, variables, mc_thresholds) {

  stat_fns <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)
  year.num <- length(unique(format(data$date, "%Y")))

  stats.season <- .summarize_grouped_stats(
    df = data,
    variables = variables,
    group_vars = c("id", "mon"),
    stat_fns = stat_fns
  )

  stats.mon.aavg <- .summarize_grouped_stats(
    df = data,
    variables = variables,
    group_vars = c("year", "mon"),
    stat_fns = stat_fns
  )

  stats.annual.aavg <- .summarize_grouped_stats(
    df = data,
    variables = variables,
    group_vars = "year",
    stat_fns = stat_fns
  ) %>%
    dplyr::mutate(year = .data$year - min(.data$year) + 1)

  wetdry <- data %>%
    dplyr::left_join(mc_thresholds, by = c("id", "mon")) %>%
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
        threshold <- if (!is.finite(th) || th <= 0) 0 else th
        spells <- compute_spell_lengths(.data$precip, threshold = threshold, below = TRUE)
        if (length(spells) == 0) NA_real_ else mean(spells)
      },
      Wet_spells = {
        th <- .data$wet.th[1]
        threshold <- if (!is.finite(th) || th <= 0) 0 else th
        spells <- compute_spell_lengths(.data$precip, threshold = threshold, below = FALSE)
        if (length(spells) == 0) NA_real_ else mean(spells)
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


  cor.data <- .compute_anomaly_correlations(data, variables)

  cor.cond <- .compute_conditional_precip_correlations(
    data = data,
    variables = variables,
    mc_thresholds = mc_thresholds,
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

#' Summarize grouped statistics
#'
#' @param df Data frame to summarize.
#' @param variables Character vector of variables to summarize.
#' @param group_vars Character vector of grouping columns.
#' @param stat_fns Named list of summary functions.
#' @return Long-format data frame with `variable`, `stat`, and `value`.
#' @keywords internal
.summarize_grouped_stats <- function(df, variables, group_vars, stat_fns) {

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarize(
      dplyr::across(dplyr::all_of(variables), stat_fns, .names = "{.col}:{.fn}"),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(group_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    tidyr::separate(.data$variable, into = c("variable", "stat"), sep = ":") %>%
    dplyr::mutate(
      stat = factor(.data$stat, levels = names(stat_fns)),
      variable = as.character(.data$variable)
    )
}

#' Compute correlation matrix on monthly anomalies
#'
#' @param data Data frame with `date`, `id`, and variables.
#' @param variables Character vector of variable names to correlate.
#' @return Long-format data frame of upper-triangle pairwise correlations.
#' @keywords internal
.compute_anomaly_correlations <- function(data, variables) {

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
#'
#' @param data Data frame with daily values and `date`/`id` columns.
#' @param variables Character vector of variable names (must include `"precip"`).
#' @param mc_thresholds Optional wet-day thresholds (required for `wet_def = "monthly_quantile"`).
#' @param wet_def Wet-day definition: `"gt0"` or `"monthly_quantile"`.
#' @param use_anom Logical; subtract monthly means before correlation.
#' @param use_log_precip_on_wet Logical; apply `log1p` to precip on wet days.
#' @param method Correlation method passed to `stats::cor`.
#' @param min_pairs Minimum number of paired values required to report a correlation.
#' @return Data frame of conditional correlations by grid, regime, and variable.
#' @keywords internal
.compute_conditional_precip_correlations <- function(data,
                                           variables,
                                           mc_thresholds = NULL,
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
    if (is.null(mc_thresholds)) stop("mc_thresholds required when wet_def='monthly_quantile'", call. = FALSE)
    dat_wide <- dat_wide %>% dplyr::left_join(mc_thresholds, by = c("id", "mon"))
    wet_flag <- dat_wide$precip >= dat_wide$wet.th
  } else {
    wet_flag <- dat_wide$precip > 0
  }
  dry_flag <- !wet_flag & !is.na(dat_wide$precip)

  # internal function for safe cor
  .safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if (n < min_pairs) return(c(value = NA_real_, n = n))
    c(value = stats::cor(x[ok], y[ok], method = method), n = n)
  }

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


#' Filter correlations by type (cross-grid or inter-variable)
#'
#' @param cor_data Data frame of correlations.
#' @param same_var Logical; require same-variable pairs.
#' @param same_id Logical; require same-grid pairs.
#' @param allowed_pairs Character vector of allowed pair keys.
#' @param pair_type `"id"` for grid pairs or `"variable"` for variable pairs.
#' @return Filtered data frame of correlations.
#' @keywords internal
.filter_correlation_pairs <- function(cor_data, same_var, same_id, allowed_pairs, pair_type) {

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

#' Build data for plotting
#'
#' @param obs_results Output from `.summarize_observed_data()`.
#' @param sim_results Output from `.summarize_simulated_data()`.
#' @param variables Character vector of variable names.
#' @return List of standardized data frames for diagnostic plots.
#' @keywords internal
.build_plot_data <- function(obs_results, sim_results, variables) {

  # Define statistical functions for reference
  stat_fns <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)

  # Merge seasonal statistics
  daily.stats.season <- sim_results$stats.season %>%
    dplyr::left_join(
      obs_results$stats.season,
      by = c("id", "mon", "variable", "stat")
    ) %>%
    dplyr::mutate(
      variable = factor(.data$variable),
      stat = factor(.data$stat, levels = names(stat_fns))
    )

  # Merge correlations (unconditional anomalies)
  stats.allcor <- sim_results$cor %>%
    dplyr::left_join(
      obs_results$cor,
      by = c("id1", "variable1", "id2", "variable2")
    )

  # Merge conditional precip correlations (within-grid; new)
  stats.allcor.cond <- sim_results$cor.cond %>%
    dplyr::left_join(
      obs_results$cor.cond,
      by = c("id1", "variable1", "id2", "variable2", "regime", "transform")
    )

  # ---------------------------------------------------------------------------
  # Robust ID + variable pairing without assuming ids are 1..n
  # ---------------------------------------------------------------------------

  # Prefer IDs from the observed long data if available; otherwise infer from correlations
  ids_from_obs <- NULL
  if (!is.null(obs_results$data) && "id" %in% names(obs_results$data)) {
    ids_from_obs <- sort(unique(obs_results$data$id))
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

  stats.crosscor <- .filter_correlation_pairs(
    cor_data = stats.allcor,
    same_var = TRUE,
    same_id = FALSE,
    allowed_pairs = id_pairs_allowed,
    pair_type = "id"
  )

  # ---------------------------------------------------------------------------
  # Inter-variable correlations: same grid, different variables
  # ---------------------------------------------------------------------------

  stats.intercor <- .filter_correlation_pairs(
    cor_data = stats.allcor,
    same_var = FALSE,
    same_id = TRUE,
    allowed_pairs = var_pairs_allowed,
    pair_type = "variable"
  )

  # ---------------------------------------------------------------------------
  # Wet/dry statistics
  # ---------------------------------------------------------------------------

  stats.wetdry <- sim_results$wetdry %>%
    dplyr::left_join(
      obs_results$wetdry,
      by = c("id", "mon", "variable", "stat", "type")
    ) %>%
    dplyr::mutate(stat = factor(.data$stat, levels = c("Dry", "Wet")))

  list(
    daily_stats_season = daily.stats.season,
    stats_mon_aavg_sim = sim_results$stats.mon.aavg,
    stats_mon_aavg_obs = obs_results$stats.mon.aavg,
    stats_annual_aavg_sim = sim_results$stats.annual.aavg,
    stats_annual_aavg_obs = obs_results$stats.annual.aavg,
    stats_crosscor = stats.crosscor,
    stats_intercor = stats.intercor,
    stats_wetdry = stats.wetdry,
    stats_precip_cor_cond = stats.allcor.cond
  )
}






