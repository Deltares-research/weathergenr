#' Evaluate Synthetic Weather Simulations Against Historical Observations
#'
#' This function compares daily synthetic weather simulations with observed data
#' across multiple grid cells. It computes summary statistics, wet/dry day counts,
#' spell lengths, and inter-site correlations, and generates diagnostic plots to
#' assess the performance of stochastic weather generators.
#'
#' @param daily.sim A list of simulated weather realizations. Each element should be a list of data frames, one per grid cell, containing daily values and a `date` column.
#' @param daily.obs A list of observed weather data frames, one per grid cell. Each should contain a `date` column and the variables specified in `variables`.
#' @param output.path File path to save generated plots. If `NULL`, plots are not saved.
#' @param variables A character vector of variable names to evaluate (e.g., `c("precip", "temp")`).
#' @param variable.labels Optional character vector of variable labels for use in plots. Defaults to `variables` if `NULL`.
#' @param variable.units Optional character vector of variable units. Defaults to empty strings.
#' @param realization.num Integer. Number of synthetic realizations in `daily.sim`.
#' @param wet.quantile Quantile threshold for wet days (default is 0.2).
#' @param extreme.quantile Quantile threshold for extremely wet days (default is 0.8).
#' @param show.title Logical. Whether to display titles in the generated plots (default is `TRUE`).
#' @param save.plots Logical. Whether to save plots to `output.path` (default is `TRUE`).
#'
#' @importFrom stats cor setNames
#' @importFrom utils combn
#'
#' @return A named list of `ggplot2` plot objects. These include:
#' \itemize{
#'   \item \code{daily_means}, \code{daily_sd}: Daily mean and standard deviation comparisons.
#'   \item \code{spell_lengths}, \code{spell_duration}: Wet/dry spell length and frequency plots.
#'   \item \code{crossgrid_cor}, \code{intergrid_cor}: Correlation plots between sites and variables.
#'   \item \code{annual_cycle}: Monthly mean climatology.
#'   \item \code{annual_mean}: Annual mean precipitation time series.
#'   \item \code{annual_pattern_<var>}: Monthly summary statistics per variable.
#' }
#'
#' @details
#' - Observed and simulated data must be structured consistently.
#' - Precipitation variable is used for wet/dry analysis and must be included in `variables`.
#' - Computations include monthly, annual, and cross-site statistics.
#' - Extensive use of `dplyr`, `tidyr`, `ggplot2`, and `e1071` packages.
#'
#' @examples
#' \dontrun{
#' plots <- evaluateWegen(
#'   daily.sim = synthetic_data,
#'   daily.obs = observed_data,
#'   output.path = "output/",
#'   variables = c("precip", "temp"),
#'   variable.labels = c("Precipitation", "Temperature"),
#'   realization.num = 3
#' )
#' plots$daily_means
#' }
#'
#' @export
evaluateWegen <- function(
    daily.sim = NULL,
    daily.obs = NULL,
    output.path = NULL,
    variables = NULL,
    variable.labels = NULL,
    variable.units = NULL,
    realization.num = NULL,
    wet.quantile = 0.2,
    extreme.quantile = 0.8,
    show.title = TRUE,
    save.plots = TRUE) {

  # Create output directory if doesn't exist
  if (!dir.exists(output.path)) {
    dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
  }

  # ----- Helper Functions -----------------------------------------------------
  gg_multipanel_export <- function(p, p.name, show.title = TRUE, save.plots = TRUE,
                                   p.title = NULL, p.subtitle = NULL, output.path) {
    if (show.title && !is.null(p.title)) {
      p <- p + labs(title = p.title, subtitle = p.subtitle)
    }
    if (save.plots) {
      ncol <- p$facet$params$ncol
      nrow <- p$facet$params$nrow
      if (is.null(ncol)) ncol <- 2
      if (is.null(nrow)) nrow <- 2
      w <- ncol * 4
      h <- nrow * 4 + 0.5
      ggsave(file.path(output.path, p.name), plot = p, width = w, height = h)
    }
    p
  }

  stat_funs <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)
  compute_grouped_stats <- function(df, variables = NULL, group_vars, stat_funs = NULL, values_to = "value") {
    stopifnot(is.data.frame(df))
    stopifnot(all(group_vars %in% names(df)))
    if (is.null(variables)) {
      variables <- setdiff(
        names(df)[vapply(df, is.numeric, logical(1))],
        group_vars
      )
    }
    if (is.null(stat_funs)) {
      stat_funs <- list(mean = mean, sd = stats::sd, skewness = e1071::skewness)
    }
    df %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(across(all_of(variables), stat_funs, .names = "{.col}:{.fn}"), .groups = "drop") %>%
      pivot_longer(cols = -all_of(group_vars), names_to = "variable", values_to = values_to) %>%
      separate(variable, into = c("variable", "stat"), sep = ":") %>%
      mutate(stat = factor(stat, levels = names(stat_funs)),
             variable = as.character(variable))
  }

  generate_symmetric_dummy_points <- function(df, facet_var = "variable",
                                              x_col = "Observed", y_col = "Simulated") {
    df %>%
      filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]])) %>%
      filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]])) %>%
      mutate(max_val = pmax(.data[[x_col]], .data[[y_col]]),
             min_val = pmin(.data[[x_col]], .data[[y_col]])) %>%
      group_by(.data[[facet_var]]) %>%
      summarize(
        minlim = min(min_val, na.rm = TRUE),
        maxlim = max(max_val, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::expand_grid(lim = c("min", "max")) %>%
      mutate(
        Observed = ifelse(lim == "min", minlim, maxlim),
        Simulated = ifelse(lim == "min", minlim, maxlim)
      ) %>%
      select(all_of(facet_var), Observed, Simulated)
  }


  # ----- Input Checks & Setup -------------------------------------------------

  # Check variables
  if (is.null(variable.labels)) variable.labels <- variables
  if (is.null(variable.units)) variable.units <- rep("", length(variables))
  options(dplyr.summarise.inform = FALSE, tidyverse.quiet = TRUE)
  nsgrids <- length(daily.sim[[1]])
  logger::log_info("Comparison accross {nsgrids} grid cells")
  stat_level <- c("mean", "sd", "skewness")
  stat_label <- c("mean", "standard dev.", "skewness")

  # Plotting variables
  pl_sub <- "Value range and median values from all simulations are shown against the observed"
  font_size <- 12; alpha_val <- 0.4
  theme_wgplots <- theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 14), plot.subtitle = element_text(size = 10))
  plot_cols <- setNames(c("blue3", "gray40"), c("Observed", "Simulated"))
  plots <- list()


  # ----- Observed Climate Stats -----------------------------------------------

  logger::log_info("Calculating historical trace statistics")

  hist_year_num <- nrow(daily.obs[[1]]) / 365
  hist_daily_tidy <- lapply(daily.obs, "[", c("date", variables)) %>%
    bind_rows(.id = "id") %>%
    mutate(
      year = as.numeric(format(date, "%Y")),
      mon = as.numeric(format(date, "%m")),
      day = as.numeric(format(date, "%d")),
      id = as.numeric(id)
    ) %>%
    select(id, year, mon, day, all_of(variables))

  hist_stats_season <- compute_grouped_stats(hist_daily_tidy, variables, c("id", "mon"), values_to = "Observed")
  hist_stats_mon_aavg <- compute_grouped_stats(hist_daily_tidy, variables, c("year", "mon"), values_to = "Observed")
  hist_stats_annual_aavg <- compute_grouped_stats(hist_daily_tidy, variables, c("year"), values_to = "Observed") %>%
    mutate(year = year - min(year) + 1)
  hist_stats_season_aavg <- compute_grouped_stats(hist_daily_tidy, variables, c("mon"), values_to = "Observed")

  mc_thresholds <- hist_daily_tidy %>%
    group_by(year, mon, day) %>%
    summarize(precip = mean(precip), .groups = "drop") %>%
    group_by(mon) %>%
    summarize(
      wet_th = stats::quantile(precip, wet.quantile, names = FALSE),
      extreme_th = stats::quantile(precip, extreme.quantile, names = FALSE),
      .groups = "drop"
    )

  dfx <- hist_daily_tidy %>%
    left_join(mc_thresholds, by = "mon") %>%
    group_by(id, mon)

  hist_wetdry_days <- dfx %>%
    summarize(
      Wet = length(which(precip >= wet_th)) / hist_year_num,
      Dry = length(which(precip < wet_th)) / hist_year_num,
      .groups = "drop"
    ) %>%
    pivot_longer(cols = Wet:Dry, names_to = "stat", values_to = "Observed") %>%
    mutate(variable = "precip", .after = mon)

  hist_wetdry_spells <- dfx %>%
    summarize(
      dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
      wet = averageSpellLength(precip, threshold = wet_th, below = FALSE),
      .groups = "drop"
    ) %>%
    group_by(id, mon) %>%
    summarize(Dry = mean(dry), Wet = mean(wet), .groups = "drop") %>%
    pivot_longer(cols = Dry:Wet, names_to = "stat", values_to = "Observed") %>%
    mutate(variable = "precip", .after = mon)

  hist_allcor_ini <- hist_daily_tidy %>%
    pivot_longer(cols = -c(id, year, mon, day), names_to = "variable", values_to = "value") %>%
    unite(id_variable, id, variable, sep = ":") %>%
    pivot_wider(names_from = id_variable, values_from = value) %>%
    select(-year, -mon, -day) %>%
    cor(use = "pairwise.complete.obs") %>%
    as_tibble()

  hist_allcor <- hist_allcor_ini %>%
    mutate(id_variable1 = colnames(hist_allcor_ini), .before = 1) %>%
    pivot_longer(cols = -id_variable1, names_to = "id_variable2", values_to = "Observed") %>%
    separate(id_variable1, c("id1", "variable1"), sep = ":") %>%
    separate(id_variable2, c("id2", "variable2"), sep = ":")

  # ----- Simulated Climate Stats ----------------------------------------------

  logger::log_info("Calculating synthetic trace statistics")

  sim_stats_season <- vector("list", realization.num)
  sim_stats_mon_aavg <- vector("list", realization.num)
  sim_stats_annual_aavg <- vector("list", realization.num)
  sim_allcor <- vector("list", realization.num)
  sim_wetdry_days <- vector("list", realization.num)
  sim_wetdry_spells <- vector("list", realization.num)
  sim_year_num <- nrow(daily.sim[[1]][[1]]) / 365

  sim_date_tbl <- daily.sim[[1]][[1]] %>%
    select(date) %>%
    mutate(
      year = as.numeric(format(date, "%Y")),
      mon = as.numeric(format(date, "%m")),
      day = as.numeric(format(date, "%d"))
    )

  # Loop over each realization and calculate statistics
  for (n in 1:realization.num) {

    sim_daily_tidy <- lapply(daily.sim[[n]], "[", c("date", variables)) %>%
      bind_rows(.id = "id") %>%
      left_join(sim_date_tbl, by = "date") %>%
      select(id, year, mon, day, all_of(variables))

    sim_stats_season[[n]] <- compute_grouped_stats(sim_daily_tidy, variables, c("id", "mon"), values_to = "Simulated")
    sim_stats_mon_aavg[[n]] <- compute_grouped_stats(sim_daily_tidy, variables, c("year", "mon"), values_to = "Simulated")
    sim_stats_annual_aavg[[n]] <- compute_grouped_stats(sim_daily_tidy, variables, c("year"), values_to = "Simulated") %>%
      mutate(year = year - min(year) + 1)

    sim_allcor_ini <- sim_daily_tidy %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value") %>%
      unite(id_variable, c("id", "variable"), sep = ":") %>%
      pivot_wider(names_from = id_variable, values_from = value) %>%
      select(-year, -mon, -day) %>%
      cor(use = "pairwise.complete.obs") %>%
      as_tibble()

    sim_allcor[[n]] <- sim_allcor_ini %>%
      mutate(id_variable1 = colnames(sim_allcor_ini), .before = 1) %>%
      pivot_longer(-id_variable1, names_to = "id_variable2", values_to = "Simulated") %>%
      separate(id_variable1, c("id1", "variable1"), sep = ":") %>%
      separate(id_variable2, c("id2", "variable2"), sep = ":")

    dfx2 <- sim_daily_tidy %>%
      left_join(mc_thresholds, by = "mon") %>%
      group_by(id, mon)

    sim_wetdry_days[[n]] <- dfx2 %>%
      summarize(
        Wet = length(which(precip >= wet_th)) / sim_year_num,
        Dry = length(which(precip < wet_th)) / sim_year_num,
        .groups = "drop"
      ) %>%
      pivot_longer(cols = Wet:Dry, names_to = "stat", values_to = "value") %>%
      mutate(variable = "precip") %>%
      select(id, mon, variable, stat, Simulated = value)

    sim_wetdry_spells[[n]] <- dfx2 %>%
      summarize(
        Dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
        Wet = averageSpellLength(precip, threshold = wet_th, below = FALSE),
        .groups = "drop"
      ) %>%
      group_by(id, mon) %>%
      summarize(Dry = mean(Dry), Wet = mean(Wet), .groups = "drop") %>%
      pivot_longer(cols = Dry:Wet, names_to = "stat", values_to = "Simulated") %>%
      mutate(variable = "precip", .before = stat)
  }

  sim_stats_season <- bind_rows(sim_stats_season, .id = "rlz") %>% mutate(id = as.numeric(id))
  sim_stats_mon_aavg <- bind_rows(sim_stats_mon_aavg, .id = "rlz")
  sim_stats_annual_aavg <- bind_rows(sim_stats_annual_aavg, .id = "rlz") %>% mutate(year = year - min(year) + 1)
  sim_allcor <- bind_rows(sim_allcor, .id = "rlz")
  sim_wetdry_days <- bind_rows(sim_wetdry_days, .id = "rlz") %>% mutate(id = as.numeric(id))
  sim_wetdry_spells <- bind_rows(sim_wetdry_spells, .id = "rlz") %>% mutate(id = as.numeric(id))

  # ----- Merge Results --------------------------------------------------------

  logger::log_info("Merging statistics from observed and simulated series")

  var_combs <- apply(combn(variables, 2), 2, paste, collapse = ":")
  id_combs <- apply(combn(1:nsgrids, 2), 2, paste, collapse = ":")
  daily_stats_season <- sim_stats_season %>%
    left_join(hist_stats_season, by = c("id", "mon", "variable", "stat")) %>%
    mutate(variable = factor(variable), stat = factor(stat, levels = names(stat_funs)))
  stats_allcor <- sim_allcor %>%
    left_join(hist_allcor, by = c("id1", "variable1", "id2", "variable2"))
  stats_intersite_cor <- stats_allcor %>%
    filter(variable1 == variable2, id1 != id2) %>%
    unite(id, c("id1", "id2"), sep = ":") %>%
    filter(id %in% id_combs)
  stats_cross_cor <- stats_allcor %>%
    filter(id1 == id2, variable1 != variable2) %>%
    unite(id, c("id1", "id2"), sep = ":") %>%
    unite(variable, c("variable1", "variable2"), sep = ":") %>%
    filter(variable %in% var_combs)
  stats_wetdry_days <- sim_wetdry_days %>%
    left_join(hist_wetdry_days, by = c("id", "mon", "variable", "stat")) %>%
    mutate(stat = factor(stat, levels = c("Dry", "Wet")))
  stats_wetdry_spells <- sim_wetdry_spells %>%
    left_join(hist_wetdry_spells, by = c("id", "mon", "variable", "stat")) %>%
    mutate(stat = factor(stat, levels = c("Dry", "Wet")))

  # ----- Plot Results ---------------------------------------------------------
  logger::log_info("Preparing comparison plots")

  # 1) Daily mean statistics for all variables
  dummy_gg <- generate_symmetric_dummy_points(
    df = filter(daily_stats_season, stat == "mean"),
    facet_var = "variable", x_col = "Observed", y_col = "Simulated")

  p <- ggplot(
    filter(daily_stats_season, stat == "mean"),
    aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    # These invisible points define symmetric axes in each panel
    geom_point(data = dummy_gg, aes(x = Observed, y = Simulated), color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange", fun.max = max, fun.min = min,
      alpha = alpha_val, linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = 2, nrow = 2)


  plots[[1]] <- gg_multipanel_export(p = p,
                                     p.name = "daily_mean.png",
                                     p.title = "Daily means for all grid cell and months",
                                     p.subtitle = pl_sub,
                                     save.plots = save.plots,
                                     output.path = output.path)


  # 2) Daily standard deviations for all variables
  dummy_gg <- generate_symmetric_dummy_points(
    df = filter(daily_stats_season, stat == "sd"),
    facet_var = "variable", x_col = "Observed", y_col = "Simulated")

  plot_name <- "daily_sd.png"

  p <- ggplot(
    filter(daily_stats_season, stat == "sd"),
    aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_point(data = dummy_gg, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange", fun.max = max, fun.min = min,
      alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = 2, nrow = 2)

  plots[[2]] <- gg_multipanel_export(p = p,
                                     p.name = "daily_sd.png",
                                     p.title = "Daily standard deviations for all grid cell and months",
                                     p.subtitle = pl_sub,
                                     save.plots = save.plots,
                                     output.path = output.path)

  # 3) Wet and dry spell statistics
  dummy_gg <- stats_wetdry_spells %>%
    group_by(stat) %>%
    summarize(
      minval = min(Simulated, Observed) * 1,
      maxval = max(Simulated, Observed) * 1) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(stat, Observed = value, Simulated = value)

  p <- ggplot(stats_wetdry_spells, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(
      geom = "linerange", fun.max = max, fun.min = min,
      alpha = alpha_val, linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(stat ~ ., ncol = 2, nrow =1 , scales = "free") +
    geom_point(data = dummy_gg, color = "blue", alpha = 0) +
    xlab("Observed") +
    ylab("Simulated")

  plots[[3]] <- gg_multipanel_export(p = p,
                                     p.name = "spell_length.png",
                                     p.title = "Average dry and wet spell length per month, across all grid cells",
                                     p.subtitle = pl_sub,
                                     save.plots = save.plots,
                                     output.path = output.path)

  # 4) Average number of wet and dry days
  dummy_gg <- stats_wetdry_days %>%
    group_by(stat) %>%
    summarize(
      minval = min(Simulated, Observed) * 1,
      maxval = max(Simulated, Observed) * 1
    ) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(stat, Observed = value, Simulated = value)

  p <- ggplot(stats_wetdry_days, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    geom_point(data = dummy_gg, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange", fun.max = max, fun.min = min,
      alpha = alpha_val, linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(stat ~ ., ncol = 2, nrow = 1, scales = "free") +
    xlab("Observed") +
    ylab("Simulated")

  plots[[4]] <- gg_multipanel_export(p = p,
                                     p.name = "wetdry_days_count.png",
                                     p.title = "Average number of dry and wet days per month accross all grid cells",
                                     p.subtitle = pl_sub,
                                     save.plots = save.plots,
                                     output.path = output.path)

  # 5) INTERGRID CORRELATIONS
  dummy_gg <- generate_symmetric_dummy_points(
    df = stats_intersite_cor,
    facet_var = "variable1", x_col = "Observed", y_col = "Simulated")

  p <- ggplot(stats_intersite_cor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    geom_point(data = dummy_gg, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange", fun.max = max, fun.min = min,
      alpha = alpha_val, linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(variable1 ~ ., ncol = 2, nrow = 2, scales = "free") +
    xlab("Observed") +
    ylab("Simulated")

  plots[[5]] <- gg_multipanel_export(p = p,
                                     p.name = "intergrid_correlations.png",
                                     p.title = "Inter-grid correlations",
                                     p.subtitle = paste0(pl_sub, "\nCorrelations are calculated over daily series"),
                                     save.plots = save.plots,
                                     output.path = output.path)



  # 6) CROSS-GRID CORRELATIONS
  dummy_gg <- generate_symmetric_dummy_points(
    df = stats_cross_cor,
    facet_var = "variable", x_col = "Observed", y_col = "Simulated")

  p <- ggplot(stats_cross_cor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    geom_point(data = dummy_gg, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange", fun.max = max, fun.min = min,
      alpha = alpha_val, linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(variable ~ ., ncol = 3, nrow = 2, scales = "free") +
    xlab("Observed") +
    ylab("Simulated")

  plots[[6]] <- gg_multipanel_export(p = p,
                                     p.name = "crossgrid_correlations.png",
                                     p.title = "Cross-grid correlations",
                                     p.subtitle = paste0(pl_sub, "\nCorrelations are calculated over daily series"),
                                     save.plots = save.plots,
                                     output.path = output.path)


  # 7) Monthly statistics per variable
  for (v in 1:length(variables)) {
    dat <- sim_stats_mon_aavg %>%
      filter(variable == variables[v]) %>%
      mutate(type = "Simulated") %>%
      rename(value = Simulated)

    dat2 <- hist_stats_mon_aavg %>%
      mutate(rlz = "0", .before = year) %>%
      filter(variable == variables[v] & !is.nan(value)) %>%
      mutate(type = "Observed") %>%
      rename(value = Observed)

    datx <- bind_rows(dat, dat2)

    p <- ggplot(datx, aes(x = as.factor(mon), y = value, fill = type, color = type)) +
      theme_wgplots +
      geom_boxplot(alpha = 0.2) +
      facet_wrap(~stat, scales = "free", ncol = 2, nrow = 2) +
      scale_fill_manual("", values = plot_cols) +
      scale_color_manual("", values = plot_cols) +
      stat_summary(fun = "mean", size = 3, geom = "point", position = position_dodge(0.8), shape = 18) +
      labs(x = "", y = "") +
      theme(
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 12)
      ) +
      scale_x_discrete(labels = substr(month.name, 1, 1))

    plots[[5+v]] <- gg_multipanel_export(p = p,
                                         p.name = paste0("annual_pattern_", variables[v], ".png"),
                                         p.title = paste0("Monthly patterns for ", variable.labels[v]),
                                         p.subtitle = paste0(pl_sub, "\nResults are averaged accross all grid cells."),
                                         save.plots = save.plots,
                                         output.path = output.path)

  }

  # 8) Monthly mean cycle
  sim_stats_season_aavg <- sim_stats_season %>%
    group_by(rlz, mon, variable) %>%
    filter(stat == "mean") %>%
    summarize(value = mean(Simulated)) %>%
    mutate(type = "Simulated")

  hist_stats_season_aavg2 <- hist_stats_season %>%
    group_by(mon, variable) %>%
    filter(stat == "mean") %>%
    summarize(value = mean(Observed)) %>%
    mutate(type = "Observed")

  p <- ggplot(sim_stats_season_aavg, aes(x = as.factor(mon), y = value)) +
    theme_wgplots +
    facet_wrap(~variable, scales = "free", ncol = 2, nrow = 2) +
    geom_line(aes(group = rlz, color = rlz), alpha = 0.8) +
    geom_line(data = hist_stats_season_aavg2, color = "black", group = 1, linewidth = 1.25) +
    scale_x_discrete(labels = substr(month.name, 1, 1)) +
    labs(x = "", y = "")

  plots[[length(plots)+1]] <- gg_multipanel_export(p = p,
                                                   p.name = "monthly_cycle.png",
                                                   p.title = paste0("Annual cycles of variables"),
                                                   p.subtitle =  paste0(pl_sub, "\nResults are averaged accross each month"),
                                                   save.plots = save.plots,
                                                   output.path = output.path)

  # Annual precip means as time-series
  sim_annual_aavg_precip <- sim_stats_annual_aavg %>%
    filter(stat == "mean") %>%
    filter(variable == "precip")
  hist_annual_aavg_precip <- hist_stats_annual_aavg %>%
    filter(stat == "mean") %>%
    filter(variable == "precip")

  plots[[length(plots)+1]] <- ggplot(sim_annual_aavg_precip, aes(x = year)) +
    theme_wgplots +
    geom_line(aes(y = Simulated, group = rlz), color = "gray30", alpha = 0.4) +
    geom_point(aes(y = Simulated, group = rlz), size = 0.5, color = "gray30", alpha = 0.3) +
    geom_line(aes(y = Observed),
              data = hist_annual_aavg_precip, color = "blue", group = 1) +
    geom_point(aes(y = Observed),
               size = 0.5, alpha = 0.3,
               data = hist_annual_aavg_precip, color = "blue", group = 1
    ) +
    labs(x = "serial year", y = "mm/day")

  if (show.title) {
    plots[[length(plots)]] <- plots[[length(plots)]] +
      labs(title = paste0("Annual mean precipitation"))
  }

  if (save.plots) {
    ggsave(file.path(output.path, paste0("annual_precip.png")), height = 4, width = 8)
  }

  logger::log_info("Evaluation completed!")

  plots
}
