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
#'
#' @return A named list of `ggplot2` plot objects with class "weather_assessment".
#'   Includes:
#'   \itemize{
#'     \item \code{daily_mean}, \code{daily_sd}: Daily statistics comparisons
#'     \item \code{spell_length}, \code{wetdry_days_count}: Spell diagnostics
#'     \item \code{crossgrid}, \code{intergrid}: Correlation diagnostics
#'     \item \code{monthly_cycle}: Monthly climatology
#'     \item \code{annual_precip}: Annual precipitation time series
#'     \item \code{annual_pattern_<var>}: Per-variable monthly statistics
#'   }
#'
#' @details
#' The function performs comprehensive validation of weather generator output by
#' comparing simulated and observed statistics across multiple dimensions:
#' temporal (daily, monthly, annual), spatial (cross-grid correlations), and
#' distributional (means, variances, spell lengths).
#'
#' Observed and simulated data must be structurally consistent. The precipitation
#' variable is required for wet/dry occurrence analysis. Computations leverage
#' `dplyr`, `tidyr`, `ggplot2`, and `e1071` packages.
#'
#' @examples
#' \dontrun{
#' # Assess weather generator with 3 realizations
#' assessment <- evaluate_weather_generator(
#'   daily.sim = synthetic_data,
#'   daily.obs = observed_data,
#'   variables = c("precip", "temp"),
#'   variable.labels = c("Precipitation", "Temperature"),
#'   realization.num = 3,
#'   output.path = "diagnostics/"
#' )
#'
#' # View specific diagnostic
#' print(assessment$daily_mean)
#' }
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
    max.grids = 25) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  if (!is.numeric(max.grids) || length(max.grids) != 1 || max.grids < 1) {
    stop("'max.grids' must be a positive integer")
  }

  validate_inputs(
    daily.sim = daily.sim,
    daily.obs = daily.obs,
    variables = variables,
    realization.num = realization.num,
    wet.quantile = wet.quantile,
    extreme.quantile = extreme.quantile
  )

  # ============================================================================
  # SETUP
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(
      paste0(
        "[Assessment] Start | grids={length(daily.obs)} | ",
        "realizations={realization.num} | ",
        "variables={paste(variables, collapse = ',')} | ",
        "wet.q={wet.quantile} | extreme.q={extreme.quantile}"
      )
    )
  }
  # Set defaults
  if (is.null(variable.labels)) variable.labels <- variables

  # Output directory handling
  if (!is.null(output.path)) {
    if (!dir.exists(output.path)) {
      dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    save.plots <- FALSE
  }

  # Suppress dplyr messages
  options(dplyr.summarise.inform = FALSE, tidyverse.quiet = TRUE)

  # Plot configuration
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

  # Results container
  plots <- list()

  # ============================================================================
  # GRID SUBSAMPLING TO CONTROL MEMORY USE
  # ============================================================================

  # Grid dimensions
  n_grids <- length(daily.obs)
  n_grids_org <- n_grids

  if (n_grids > max.grids) {

    sel.grids <- sort(sample(seq_len(n_grids), max.grids))

    # Subset observed data
    daily.obs <- daily.obs[sel.grids]

    # Subset simulated data (each realization)
    daily.sim <- lapply(daily.sim, function(rlz) {
      rlz[sel.grids]
    })

    n_grids <- length(daily.obs)

    if (requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn(
        "[Assessment] Grid count reduced from {n_grids_org} to {n_grids} for memory control"
      )
    }
  }

  # Logging
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Evaluating {n_grids} grid cells with {realization.num} realizations")
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
    logger::log_info("[Assessment] Processing simulated data")
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
    logger::log_info("[Assessment] Preparing diagnostic data")
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

  # ============================================================================
  # FINALIZE AND RETURN
  # ============================================================================

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[Assessment] Completed. {length(plots)} diagnostic plots generated")
    if (save.plots) {
      logger::log_info("[Assessment] Plots saved to: {output.path}")
    }
  }

  # Add metadata and class
  structure(
    plots,
    class = c("weather_assessment", "list"),
    metadata = list(
      n_grids = n_grids,
      realization.num = realization.num,
      variables = variables,
      assessment.date = Sys.Date()
    )
  )
}


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Validate inputs for weather assessment
#' @keywords internal
validate_inputs <- function(daily.sim, daily.obs, variables, realization.num,
                           wet.quantile, extreme.quantile) {

  # Check nulls
  if (is.null(daily.sim)) stop("'daily.sim' must not be NULL")
  if (is.null(daily.obs)) stop("'daily.obs' must not be NULL")
  if (is.null(variables)) stop("'variables' must not be NULL")
  if (is.null(realization.num)) stop("'realization.num' must not be NULL")

  # Check types
  if (!is.list(daily.sim)) stop("'daily.sim' must be a list")
  if (!is.list(daily.obs)) stop("'daily.obs' must be a list")
  if (!is.character(variables)) stop("'variables' must be a character vector")
  if (!is.numeric(realization.num) || realization.num < 1) {
    stop("'realization.num' must be a positive integer")
  }

  # Check dimensions
  if (length(daily.sim) != realization.num) {
    stop("Length of 'daily.sim' must equal 'realization.num'")
  }

  if (length(daily.obs) == 0) {
    stop("'daily.obs' must contain at least one grid cell")
  }

  # Check for required precipitation variable
  if (!"precip" %in% variables) {
    stop("'variables' must include 'precip' for wet/dry spell analysis")
  }

  # Check variables exist in data
  missing.vars <- setdiff(variables, names(daily.obs[[1]]))
  if (length(missing.vars) > 0) {
    stop("Variables not found in daily.obs: ", paste(missing.vars, collapse = ", "))
  }

  # Check date column
  if (!"date" %in% names(daily.obs[[1]])) {
    stop("'daily.obs' must contain a 'date' column")
  }

  if (!"date" %in% names(daily.sim[[1]][[1]])) {
    stop("'daily.sim' must contain a 'date' column")
  }

  # Check quantiles
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

  # Extract dates and handle leap days
  his.date <- daily.obs[[1]]$date
  leap.idx <- find_leap_days(his.date)

  if (!is.null(leap.idx)) {
    his.date <- his.date[-leap.idx]
  }

  # Create date matrix
  his.datemat <- dplyr::tibble(
    date = his.date,
    year = as.integer(format(date, "%Y")),
    mon = as.integer(format(date, "%m")),
    day = as.integer(format(date, "%d"))
  )

  # Combine all grid cells
  his <- lapply(seq_len(n_grids), function(i) {
    df <- daily.obs[[i]][, variables, drop = FALSE]
    if (!is.null(leap.idx)) {
      df <- df[-leap.idx, , drop = FALSE]
    }
    dplyr::bind_cols(his.datemat, df)
  }) %>%
    dplyr::bind_rows(.id = "id") %>%
    dplyr::mutate(id = as.integer(.data$id))

  # Compute wet/extreme thresholds
  mc.thresholds <- his %>%
    dplyr::group_by(.data$id, .data$mon) %>%
    dplyr::summarize(
      wet.th = stats::quantile(.data$precip, wet.quantile, names = FALSE, na.rm = TRUE),
      extreme.th = stats::quantile(.data$precip, extreme.quantile, names = FALSE, na.rm = TRUE),
      .groups = "drop"
    )

  # Compute statistics
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
    cor = his.stats$cor %>% dplyr::rename(Observed = .data$value)
  )
}


#' Process simulated weather data
#' @keywords internal
process_simulated_data <- function(daily.sim, realization.num, variables, mc.thresholds) {

  # Create date matrix for simulations
  sim.datemat <- dplyr::tibble(
    date = daily.sim[[1]][[1]]$date,
    year = as.integer(format(date, "%Y")),
    mon = as.integer(format(date, "%m")),
    day = as.integer(format(date, "%d"))
  )

  # Combine simulations
  sim <- lapply(seq_len(realization.num), function(i) {
    daily.sim[[i]] %>%
      dplyr::bind_rows(.id = "id") %>%
      dplyr::mutate(id = as.integer(.data$id)) %>%
      dplyr::left_join(sim.datemat, by = "date")
  })

  # Compute statistics for each realization
  sim.stats.list <- lapply(seq_len(realization.num), function(i) {
    compute_timeseries_statistics(
      data = sim[[i]],
      variables = variables,
      mc.thresholds = mc.thresholds
    )
  })

  # Combine results
  list(
    stats.season = dplyr::bind_rows(
      lapply(sim.stats.list, `[[`, "stats.season"),
      .id = "rlz"
    ) %>%
      dplyr::mutate(id = as.numeric(.data$id)) %>%
      dplyr::rename(Simulated = .data$value),

    stats.mon.aavg = dplyr::bind_rows(
      lapply(sim.stats.list, `[[`, "stats.mon.aavg"),
      .id = "rlz"
    ) %>%
      dplyr::rename(Simulated = .data$value),

    stats.annual.aavg = dplyr::bind_rows(
      lapply(sim.stats.list, `[[`, "stats.annual.aavg"),
      .id = "rlz"
    ) %>%
      dplyr::mutate(year = .data$year - min(.data$year) + 1) %>%
      dplyr::rename(Simulated = .data$value),

    cor = dplyr::bind_rows(
      lapply(sim.stats.list, `[[`, "cor"),
      .id = "rlz"
    ) %>%
      dplyr::rename(Simulated = .data$value),

    wetdry = dplyr::bind_rows(
      lapply(sim.stats.list, `[[`, "wetdry"),
      .id = "rlz"
    ) %>%
      dplyr::mutate(id = as.numeric(.data$id)) %>%
      dplyr::rename(Simulated = .data$value)
  )
}


#' Compute time series statistics
#' @keywords internal
compute_timeseries_statistics <- function(data, variables, mc.thresholds) {

  # Define statistical functions
  stat.funs <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)

  # Number of years
  year.num <- length(unique(format(data$date, "%Y")))

  # Seasonal statistics (by id and month)
  stats.season <- compute_grouped_statistics(
    df = data,
    variables = variables,
    group.vars = c("id", "mon"),
    stat.funs = stat.funs
  )

  # Monthly area-averaged statistics
  stats.mon.aavg <- compute_grouped_statistics(
    df = data,
    variables = variables,
    group.vars = c("year", "mon"),
    stat.funs = stat.funs
  )

  # Annual area-averaged statistics
  stats.annual.aavg <- compute_grouped_statistics(
    df = data,
    variables = variables,
    group.vars = "year",
    stat.funs = stat.funs
  ) %>%
    dplyr::mutate(year = .data$year - min(.data$year) + 1)

  # Wet/dry spell statistics
  wetdry <- data %>%
    dplyr::left_join(mc.thresholds, by = c("id", "mon")) %>%
    dplyr::group_by(.data$id, .data$mon) %>%
    dplyr::summarize(
      Wet_days = sum(.data$precip >= .data$wet.th) / year.num,
      Dry_days = sum(.data$precip < .data$wet.th) / year.num,
      Dry_spells = mean_spell_length(.data$precip, threshold = .data$wet.th[1], below = TRUE),
      Wet_spells = mean_spell_length(.data$precip, threshold = .data$wet.th[1], below = FALSE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(.data$Wet_days, .data$Dry_days, .data$Dry_spells, .data$Wet_spells),
      names_to = "stat.full",
      values_to = "value"
    ) %>%
    tidyr::separate(.data$stat.full, into = c("stat", "type"), sep = "_") %>%
    dplyr::mutate(variable = "precip", .after = .data$mon)

  # Correlation matrix
  cor.data <- compute_correlation_matrix(data, variables)

  list(
    stats.season = stats.season,
    stats.mon.aavg = stats.mon.aavg,
    stats.annual.aavg = stats.annual.aavg,
    wetdry = wetdry,
    cor = cor.data
  )
}


#' Compute grouped statistics
#' @keywords internal
compute_grouped_statistics <- function(df, variables, group.vars, stat.funs) {

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group.vars))) %>%
    dplyr::summarize(
      dplyr::across(
        dplyr::all_of(variables),
        stat.funs,
        .names = "{.col}:{.fn}"
      ),
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


#' Compute correlation matrix efficiently
#' @keywords internal
compute_correlation_matrix <- function(data, variables) {

  # Create wide format for correlation
  mat <- data %>%
    dplyr::select(-.data$year, -.data$mon, -.data$day) %>%
    tidyr::pivot_longer(
      cols = -c(.data$id, .data$date),
      names_to = "variable",
      values_to = "value"
    ) %>%
    tidyr::unite("id.variable", .data$id, .data$variable, sep = ":") %>%
    tidyr::pivot_wider(
      names_from = .data$id.variable,
      values_from = .data$value
    ) %>%
    dplyr::select(-.data$date) %>%
    as.matrix()

  # Check for sufficient non-NA data
  if (all(is.na(mat))) {
    warning("Correlation matrix contains only NA values")
    return(dplyr::tibble(
      id1 = character(0),
      variable1 = character(0),
      id2 = character(0),
      variable2 = character(0),
      value = numeric(0)
    ))
  }

  # Compute correlation
  cmat <- stats::cor(mat, use = "pairwise.complete.obs")

  # Extract upper triangle
  tri <- upper.tri(cmat, diag = FALSE)
  idx.i <- row(cmat)[tri]
  idx.j <- col(cmat)[tri]

  # Build correlation table
  dplyr::tibble(
    id.variable1 = colnames(cmat)[idx.i],
    id.variable2 = colnames(cmat)[idx.j],
    value = cmat[tri]
  ) %>%
    tidyr::separate(.data$id.variable1, c("id1", "variable1"), sep = ":") %>%
    tidyr::separate(.data$id.variable2, c("id2", "variable2"), sep = ":") %>%
    dplyr::arrange(.data$id1, .data$variable1, .data$id2, .data$variable2)
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

  # Merge correlations
  stats.allcor <- sim.results$cor %>%
    dplyr::left_join(
      obs.results$cor,
      by = c("id1", "variable1", "id2", "variable2")
    )

  # Cross-grid correlations (same variable, different grids)
  var.combs <- apply(utils::combn(variables, 2), 2, paste, collapse = ":")
  id.combs <- apply(utils::combn(seq_len(n_grids), 2), 2, paste, collapse = ":")

  stats.crosscor <- stats.allcor %>%
    dplyr::filter(.data$variable1 == .data$variable2, .data$id1 != .data$id2) %>%
    tidyr::unite("id", c(.data$id1, .data$id2), sep = ":") %>%
    dplyr::filter(.data$id %in% id.combs)

  # Inter-grid correlations (same grid, different variables)
  stats.intercor <- stats.allcor %>%
    dplyr::filter(.data$id1 == .data$id2, .data$variable1 != .data$variable2) %>%
    tidyr::unite("id", c(.data$id1, .data$id2), sep = ":") %>%
    tidyr::unite("variable", c(.data$variable1, .data$variable2), sep = ":") %>%
    dplyr::filter(.data$variable %in% var.combs)

  # Wet/dry statistics
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
    stats.wetdry = stats.wetdry
  )
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


#' Create all diagnostic plots
#' @keywords internal
create_all_diagnostic_plots <- function(plot.data, plot.config, variables,
                                       variable.labels, show.title,
                                       save.plots, output.path) {

  plots <- list()

  # 1. Daily means
  plots$daily_mean <- create_daily_mean_plot(
    plot.data$daily.stats.season,
    plot.config,
    show.title,
    save.plots,
    output.path
  )

  # 2. Daily standard deviations
  plots$daily_sd <- create_daily_sd_plot(
    plot.data$daily.stats.season,
    plot.config,
    show.title,
    save.plots,
    output.path
  )

  # 3. Spell lengths
  plots$spell_length <- create_spell_length_plot(
    plot.data$stats.wetdry,
    plot.config,
    show.title,
    save.plots,
    output.path
  )

  # 4. Wet/dry day counts
  plots$wetdry_days_count <- create_wetdry_days_plot(
    plot.data$stats.wetdry,
    plot.config,
    show.title,
    save.plots,
    output.path
  )

  # 5. Cross-grid correlations
  plots$crossgrid <- create_crossgrid_cor_plot(
    plot.data$stats.crosscor,
    plot.config,
    show.title,
    save.plots,
    output.path
  )

  # 6. Inter-grid correlations
  plots$intergrid <- create_intergrid_cor_plot(
    plot.data$stats.intercor,
    plot.config,
    show.title,
    save.plots,
    output.path
  )

  # 7. Monthly patterns per variable
  for (i in seq_along(variables)) {
    plot.name <- paste0("annual_pattern_", variables[i])
    plots[[plot.name]] <- create_monthly_pattern_plot(
      stats.mon.aavg.sim = plot.data$stats.mon.aavg.sim,
      stats.mon.aavg.obs = plot.data$stats.mon.aavg.obs,
      variable = variables[i],
      variable.label = variable.labels[i],
      plot.config = plot.config,
      show.title = show.title,
      save.plots = save.plots,
      output.path = output.path
    )
  }

  # 8. Monthly cycle
  plots$monthly_cycle <- create_monthly_cycle_plot(
    plot.data$daily.stats.season,
    plot.config,
    show.title,
    save.plots,
    output.path
  )

  # 9. Annual precipitation
  plots$annual_precip <- create_annual_precip_plot(
    stats.annual.aavg.sim = plot.data$stats.annual.aavg.sim,
    stats.annual.aavg.obs = plot.data$stats.annual.aavg.obs,
    plot.config = plot.config,
    show.title = show.title,
    save.plots = save.plots,
    output.path = output.path
  )

  plots
}


#' Export multi-panel plot
#' @keywords internal
export_multipanel_plot <- function(p, filename, show.title, save.plots,
                                  title = NULL, subtitle = NULL, output.path) {

  if (show.title && !is.null(title)) {
    p <- p + ggplot2::labs(title = title, subtitle = subtitle)
  }

  if (save.plots && !is.null(output.path)) {
    ncol <- p$facet$params$ncol
    nrow <- p$facet$params$nrow
    if (is.null(ncol)) ncol <- 2
    if (is.null(nrow)) nrow <- 2

    width <- ncol * 4
    height <- nrow * 4 + 0.5

    ggplot2::ggsave(
      filename = file.path(output.path, filename),
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
  }

  invisible(p)
}


#' Create daily mean plot
#' @keywords internal
create_daily_mean_plot <- function(daily.stats.season, plot.config,
                                  show.title, save.plots, output.path) {

  data.mean <- daily.stats.season %>%
    dplyr::filter(.data$stat == "mean")

  dummy.points <- generate_symmetric_dummy_points(
    df = data.mean,
    facet.var = "variable",
    x.col = "Observed",
    y.col = "Simulated"
  )

  p <- ggplot2::ggplot(data.mean, ggplot2::aes(x = .data$Observed, y = .data$Simulated)) +
    plot.config$theme +
    ggplot2::geom_point(
      data = dummy.points,
      color = "blue",
      alpha = 0
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot.config$alpha,
      linewidth = 1.5
    ) +
    ggplot2::stat_summary(
      fun = "median",
      geom = "point",
      alpha = plot.config$alpha,
      size = 2
    ) +
    ggplot2::geom_abline(color = "blue") +
    ggplot2::labs(x = "Observed", y = "Simulated") +
    ggplot2::facet_wrap(~ variable, scales = "free", ncol = 2, nrow = 2)

  export_multipanel_plot(
    p = p,
    filename = "daily_mean.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Daily means for all grid cells and months",
    subtitle = plot.config$subtitle,
    output.path = output.path
  )
}


#' Create daily standard deviation plot
#' @keywords internal
create_daily_sd_plot <- function(daily.stats.season, plot.config,
                                show.title, save.plots, output.path) {

  data.sd <- daily.stats.season %>%
    dplyr::filter(.data$stat == "sd")

  dummy.points <- generate_symmetric_dummy_points(
    df = data.sd,
    facet.var = "variable",
    x.col = "Observed",
    y.col = "Simulated"
  )

  p <- ggplot2::ggplot(data.sd, ggplot2::aes(x = .data$Observed, y = .data$Simulated)) +
    plot.config$theme +
    ggplot2::geom_point(
      data = dummy.points,
      color = "blue",
      alpha = 0
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot.config$alpha,
      linewidth = 1.5
    ) +
    ggplot2::stat_summary(
      fun = "median",
      geom = "point",
      alpha = plot.config$alpha,
      size = 2
    ) +
    ggplot2::geom_abline(color = "blue") +
    ggplot2::labs(x = "Observed", y = "Simulated") +
    ggplot2::facet_wrap(~ variable, scales = "free", ncol = 2, nrow = 2)

  export_multipanel_plot(
    p = p,
    filename = "daily_sd.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Daily standard deviations for all grid cells and months",
    subtitle = plot.config$subtitle,
    output.path = output.path
  )
}


#' Create spell length plot
#' @keywords internal
create_spell_length_plot <- function(stats.wetdry, plot.config,
                                    show.title, save.plots, output.path) {

  data.spells <- stats.wetdry %>%
    dplyr::filter(.data$type == "spells")

  dummy.points <- data.spells %>%
    dplyr::group_by(.data$stat) %>%
    dplyr::summarize(
      minval = min(.data$Simulated, .data$Observed, na.rm = TRUE),
      maxval = max(.data$Simulated, .data$Observed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(.data$minval, .data$maxval),
      names_to = "type",
      values_to = "value"
    ) %>%
    dplyr::select(.data$stat, Observed = .data$value, Simulated = .data$value)

  p <- ggplot2::ggplot(data.spells, ggplot2::aes(x = .data$Observed, y = .data$Simulated)) +
    plot.config$theme +
    ggplot2::geom_abline(color = "blue") +
    ggplot2::geom_point(
      data = dummy.points,
      color = "blue",
      alpha = 0
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot.config$alpha,
      linewidth = 1.5
    ) +
    ggplot2::stat_summary(
      fun = "median",
      geom = "point",
      alpha = plot.config$alpha,
      size = 2
    ) +
    ggplot2::facet_wrap(~ stat, ncol = 2, nrow = 1, scales = "free") +
    ggplot2::labs(x = "Observed", y = "Simulated")

  export_multipanel_plot(
    p = p,
    filename = "spell_length.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Average dry and wet spell length per month, across all grid cells",
    subtitle = plot.config$subtitle,
    output.path = output.path
  )
}


#' Create wet/dry days count plot
#' @keywords internal
create_wetdry_days_plot <- function(stats.wetdry, plot.config,
                                   show.title, save.plots, output.path) {

  data.days <- stats.wetdry %>%
    dplyr::filter(.data$type == "days")

  dummy.points <- data.days %>%
    dplyr::group_by(.data$stat) %>%
    dplyr::summarize(
      minval = min(.data$Simulated, .data$Observed, na.rm = TRUE),
      maxval = max(.data$Simulated, .data$Observed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(.data$minval, .data$maxval),
      names_to = "type",
      values_to = "value"
    ) %>%
    dplyr::select(.data$stat, Observed = .data$value, Simulated = .data$value)

  p <- ggplot2::ggplot(data.days, ggplot2::aes(x = .data$Observed, y = .data$Simulated)) +
    plot.config$theme +
    ggplot2::geom_abline(color = "blue") +
    ggplot2::geom_point(
      data = dummy.points,
      color = "blue",
      alpha = 0
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot.config$alpha,
      linewidth = 1.5
    ) +
    ggplot2::stat_summary(
      fun = "median",
      geom = "point",
      alpha = plot.config$alpha,
      size = 2
    ) +
    ggplot2::facet_wrap(~ stat, ncol = 2, nrow = 1, scales = "free") +
    ggplot2::labs(x = "Observed", y = "Simulated")

  export_multipanel_plot(
    p = p,
    filename = "wetdry_days_count.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Average number of dry and wet days per month across all grid cells",
    subtitle = plot.config$subtitle,
    output.path = output.path
  )
}


#' Create cross-grid correlation plot
#' @keywords internal
create_crossgrid_cor_plot <- function(stats.crosscor, plot.config,
                                     show.title, save.plots, output.path) {

  dummy.points <- generate_symmetric_dummy_points(
    df = stats.crosscor,
    facet.var = "variable1",
    x.col = "Observed",
    y.col = "Simulated"
  )

  p <- ggplot2::ggplot(stats.crosscor, ggplot2::aes(x = .data$Observed, y = .data$Simulated)) +
    plot.config$theme +
    ggplot2::geom_abline(color = "blue") +
    ggplot2::geom_point(
      data = dummy.points,
      color = "blue",
      alpha = 0
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot.config$alpha,
      linewidth = 1.5
    ) +
    ggplot2::stat_summary(
      fun = "median",
      geom = "point",
      alpha = plot.config$alpha,
      size = 2
    ) +
    ggplot2::facet_wrap(~ variable1, ncol = 2, nrow = 2, scales = "free") +
    ggplot2::labs(x = "Observed", y = "Simulated")

  export_multipanel_plot(
    p = p,
    filename = "crossgrid_correlations.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Cross-grid correlations",
    subtitle = paste0(plot.config$subtitle, "\nCorrelations calculated over daily series"),
    output.path = output.path
  )
}


#' Create inter-grid correlation plot
#' @keywords internal
create_intergrid_cor_plot <- function(stats.intercor, plot.config,
                                     show.title, save.plots, output.path) {

  dummy.points <- generate_symmetric_dummy_points(
    df = stats.intercor,
    facet.var = "variable",
    x.col = "Observed",
    y.col = "Simulated"
  )

  p <- ggplot2::ggplot(stats.intercor, ggplot2::aes(x = .data$Observed, y = .data$Simulated)) +
    plot.config$theme +
    ggplot2::geom_abline(color = "blue") +
    ggplot2::geom_point(
      data = dummy.points,
      color = "blue",
      alpha = 0
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot.config$alpha,
      linewidth = 1.5
    ) +
    ggplot2::stat_summary(
      fun = "median",
      geom = "point",
      alpha = plot.config$alpha,
      size = 2
    ) +
    ggplot2::facet_wrap(~ variable, ncol = 3, nrow = 2, scales = "free") +
    ggplot2::labs(x = "Observed", y = "Simulated")

  export_multipanel_plot(
    p = p,
    filename = "intergrid_correlations.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Inter-variable correlations",
    subtitle = paste0(plot.config$subtitle, "\nCorrelations calculated over daily series"),
    output.path = output.path
  )
}


#' Create monthly pattern plot for a single variable
#' @keywords internal
create_monthly_pattern_plot <- function(stats.mon.aavg.sim, stats.mon.aavg.obs,
                                       variable, variable.label, plot.config,
                                       show.title, save.plots, output.path) {

  # Simulated data
  dat.sim <- stats.mon.aavg.sim %>%
    dplyr::filter(.data$variable == !!variable) %>%
    dplyr::mutate(type = "Simulated") %>%
    dplyr::rename(value = .data$Simulated)

  # Observed data
  dat.obs <- stats.mon.aavg.obs %>%
    dplyr::mutate(rlz = "0", .before = .data$year) %>%
    dplyr::filter(.data$variable == !!variable, !is.nan(.data$Observed)) %>%
    dplyr::mutate(type = "Observed") %>%
    dplyr::rename(value = .data$Observed)

  # Combine
  dat <- dplyr::bind_rows(dat.sim, dat.obs)

  p <- ggplot2::ggplot(dat, ggplot2::aes(
    x = as.factor(.data$mon),
    y = .data$value,
    fill = .data$type,
    color = .data$type
  )) +
    plot.config$theme +
    ggplot2::geom_boxplot(alpha = 0.2) +
    ggplot2::facet_wrap(~ stat, scales = "free", ncol = 2, nrow = 2) +
    ggplot2::scale_fill_manual("", values = plot.config$colors) +
    ggplot2::scale_color_manual("", values = plot.config$colors) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 3,
      geom = "point",
      position = ggplot2::position_dodge(0.8),
      shape = 18
    ) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_x_discrete(labels = substr(month.name, 1, 1))

  export_multipanel_plot(
    p = p,
    filename = paste0("annual_pattern_", variable, ".png"),
    show.title = show.title,
    save.plots = save.plots,
    title = paste0("Monthly patterns for ", variable.label),
    subtitle = paste0(plot.config$subtitle, "\nResults averaged across all grid cells"),
    output.path = output.path
  )
}


#' Create monthly cycle plot
#' @keywords internal
create_monthly_cycle_plot <- function(daily.stats.season, plot.config,
                                     show.title, save.plots, output.path) {

  # Simulated average
  sim.avg <- daily.stats.season %>%
    dplyr::group_by(.data$rlz, .data$mon, .data$variable) %>%
    dplyr::filter(.data$stat == "mean") %>%
    dplyr::summarize(value = mean(.data$Simulated, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(type = "Simulated")

  # Observed average
  obs.avg <- daily.stats.season %>%
    dplyr::group_by(.data$mon, .data$variable) %>%
    dplyr::filter(.data$stat == "mean") %>%
    dplyr::summarize(value = mean(.data$Observed, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(type = "Observed")

  p <- ggplot2::ggplot(sim.avg, ggplot2::aes(x = as.factor(.data$mon), y = .data$value)) +
    plot.config$theme +
    ggplot2::facet_wrap(~ variable, scales = "free", ncol = 2, nrow = 2) +
    ggplot2::geom_line(ggplot2::aes(group = .data$rlz, color = .data$rlz), alpha = 0.8) +
    ggplot2::geom_line(
      data = obs.avg,
      color = "black",
      group = 1,
      linewidth = 1.25
    ) +
    ggplot2::scale_x_discrete(labels = substr(month.name, 1, 1)) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::guides(color = "none")

  export_multipanel_plot(
    p = p,
    filename = "monthly_cycle.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Annual cycles of variables",
    subtitle = paste0(plot.config$subtitle, "\nResults averaged across each month"),
    output.path = output.path
  )
}


#' Create annual precipitation plot
#' @keywords internal
create_annual_precip_plot <- function(stats.annual.aavg.sim, stats.annual.aavg.obs,
                                     plot.config, show.title, save.plots, output.path) {

  # Simulated data
  sim.precip <- stats.annual.aavg.sim %>%
    dplyr::filter(.data$stat == "mean", .data$variable == "precip")

  # Observed data
  obs.precip <- stats.annual.aavg.obs %>%
    dplyr::filter(.data$stat == "mean", .data$variable == "precip")

  p <- ggplot2::ggplot(sim.precip, ggplot2::aes(x = .data$year)) +
    plot.config$theme +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$Simulated, group = .data$rlz),
      color = "gray30",
      alpha = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$Simulated, group = .data$rlz),
      size = 0.5,
      color = "gray30",
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = obs.precip,
      ggplot2::aes(y = .data$Observed),
      color = "blue",
      group = 1
    ) +
    ggplot2::geom_point(
      data = obs.precip,
      ggplot2::aes(y = .data$Observed),
      size = 0.5,
      alpha = 0.3,
      color = "blue",
      group = 1
    ) +
    ggplot2::labs(x = "Serial year", y = "mm/day")

  if (show.title) {
    p <- p + ggplot2::labs(title = "Annual mean precipitation")
  }

  if (save.plots && !is.null(output.path)) {
    ggplot2::ggsave(
      filename = file.path(output.path, "annual_precip.png"),
      plot = p,
      height = 4,
      width = 8,
      dpi = 300
    )
  }

  invisible(p)
}
