#' Create all diagnostic plots
#'
#' @description
#' Runs the full set of diagnostic plotting routines and returns a named list of ggplot
#' objects. Optionally saves plots to disk (delegated to the individual plot exporters).
#'
#' @details
#' This helper expects the precomputed plot data returned by the evaluation pipeline.
#' It does not validate plot input structure beyond basic use in downstream plotting.
#'
#' @param plot_data List of precomputed diagnostic datasets produced by the evaluation pipeline.
#' @param plot_config List of plotting configuration options (theme, alpha, colors, subtitle).
#' @param variables Character vector of variable names to loop over for monthly pattern plots.
#' @param show_title Logical; if \code{TRUE}, titles/subtitles are added to plots where supported.
#' @param save_plots Logical; if \code{TRUE}, plots are written to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return Named list of ggplot objects for all diagnostics created.
#'
#' @examples
#' \dontrun{
#'   plot_data <- list()
#'   plot_config <- list(
#'     subtitle = "Example",
#'     alpha = 0.4,
#'     colors = c(Observed = "blue3", Simulated = "gray40"),
#'     theme = ggplot2::theme_bw()
#'   )
#'   plots <- create_all_diagnostic_plots(
#'     plot_data = plot_data,
#'     plot_config = plot_config,
#'     variables = c("precip", "temp"),
#'     show_title = FALSE,
#'     save_plots = FALSE,
#'     output_dir = NULL
#'   )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
create_all_diagnostic_plots <- function(plot_data, plot_config, variables,
                                        show_title, save_plots, output_dir) {

  plots <- list()

  plots$daily_mean <- .create_daily_mean_plot(
    daily_stats_season = plot_data$daily_stats_season,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots$daily_sd <- .create_daily_sd_plot(
    daily_stats_season = plot_data$daily_stats_season,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots$spell_length <- .create_spell_length_plot(
    stats_wetdry = plot_data$stats_wetdry,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots$wetdry_days_count <- .create_wetdry_days_plot(
    stats_wetdry = plot_data$stats_wetdry,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots$crossgrid <- .create_crossgrid_cor_plot(
    stats_crosscor = plot_data$stats_crosscor,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots$intergrid <- .create_intergrid_cor_plot(
    stats_intercor = plot_data$stats_intercor,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots$precip_cond_cor <- .create_precip_cond_cor_plot(
    stats_precip_cor_cond = plot_data$stats_precip_cor_cond,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  # Monthly patterns per variable (use var name as label)
  for (v in variables) {
    plot_name <- paste0("annual_pattern_", v)
    plots[[plot_name]] <- .create_monthly_pattern_plot(
      stats_mon_aavg_sim = plot_data$stats_mon_aavg_sim,
      stats_mon_aavg_obs = plot_data$stats_mon_aavg_obs,
      variable = v,
      plot_config = plot_config,
      show_title = show_title,
      save_plots = save_plots,
      output_dir = output_dir
    )
  }

  plots$monthly_cycle <- .create_monthly_cycle_plot(
    daily_stats_season = plot_data$daily_stats_season,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots$annual_precip <- .create_annual_precip_plot(
    stats_annual_aavg_sim = plot_data$stats_annual_aavg_sim,
    stats_annual_aavg_obs = plot_data$stats_annual_aavg_obs,
    plot_config = plot_config,
    show_title = show_title,
    save_plots = save_plots,
    output_dir = output_dir
  )

  plots
}


#' Export a faceted (multi-panel) ggplot
#'
#' Adds title/subtitle (optional) and saves the plot to disk (optional), using
#' facet layout to infer width/height. Returns the plot invisibly.
#'
#' @param p ggplot object, typically faceted via \code{facet_wrap()} or \code{facet_grid()}.
#' @param filename Character; output filename (e.g., \code{"daily_mean.png"}).
#' @param show_title Logical; if \code{TRUE}, adds \code{title}/\code{subtitle} via \code{labs()}.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param title Character; plot title (only used when \code{show_title = TRUE}).
#' @param subtitle Character; plot subtitle (only used when \code{show_title = TRUE}).
#' @param output_dir Character; output directory for saved plots.
#'
#' @return The ggplot object \code{p}, returned invisibly.
#'
#' @keywords internal
#' @import ggplot2
.export_multipanel_plot <- function(p, filename, show_title, save_plots,
                                   title = NULL, subtitle = NULL, output_dir) {

  if (show_title && !is.null(title)) {
    p <- p + labs(title = title, subtitle = subtitle)
  }

  if (save_plots && !is.null(output_dir)) {
    ncol <- p$facet$params$ncol
    nrow <- p$facet$params$nrow
    if (is.null(ncol)) ncol <- 2
    if (is.null(nrow)) nrow <- 2

    width <- ncol * 4
    height <- nrow * 4 + 0.5

    ggplot2::ggsave(
      filename = file.path(output_dir, filename),
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
  }

  invisible(p)
}


#' Create daily mean diagnostic plot
#'
#' Faceted observed-vs-simulated comparison of daily mean values by variable, using
#' summary ranges and medians across grid cells/months. Optionally saves the plot.
#'
#' @param daily_stats_season Data frame of seasonal daily statistics including columns
#'   \code{stat}, \code{Observed}, \code{Simulated}, and \code{variable}.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
.create_daily_mean_plot <- function(daily_stats_season, plot_config,
                                   show_title, save_plots, output_dir) {

  data_mean <- daily_stats_season %>%
    dplyr::filter(.data$stat == "mean")

  dummy_points <- generate_symmetric_dummy_points(
    df = data_mean,
    facet_var = "variable",
    x_col = "Observed",
    y_col = "Simulated"
  )

  p <- ggplot(data_mean, aes(x = .data$Observed, y = .data$Simulated)) +
    plot_config$theme +
    geom_point(data = dummy_points, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot_config$alpha,
      linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = plot_config$alpha, size = 2) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(~ variable, scales = "free", ncol = 2, nrow = 2)

  .export_multipanel_plot(
    p = p,
    filename = "daily_mean.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Daily means for all grid cells and months",
    subtitle = plot_config$subtitle,
    output_dir = output_dir
  )
}


#' Create daily standard deviation diagnostic plot
#'
#' Faceted observed-vs-simulated comparison of daily standard deviations by variable,
#' using summary ranges and medians across grid cells/months. Optionally saves the plot.
#'
#' @param daily_stats_season Data frame of seasonal daily statistics including columns
#'   \code{stat}, \code{Observed}, \code{Simulated}, and \code{variable}.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
.create_daily_sd_plot <- function(daily_stats_season, plot_config,
                                 show_title, save_plots, output_dir) {

  data_sd <- daily_stats_season %>%
    dplyr::filter(.data$stat == "sd")

  dummy_points <- generate_symmetric_dummy_points(
    df = data_sd,
    facet_var = "variable",
    x_col = "Observed",
    y_col = "Simulated"
  )

  p <- ggplot(data_sd, aes(x = .data$Observed, y = .data$Simulated)) +
    plot_config$theme +
    geom_point(data = dummy_points, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot_config$alpha,
      linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = plot_config$alpha, size = 2) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(~ variable, scales = "free", ncol = 2, nrow = 2)

  .export_multipanel_plot(
    p = p,
    filename = "daily_sd.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Daily standard deviations for all grid cells and months",
    subtitle = plot_config$subtitle,
    output_dir = output_dir
  )
}


#' Create wet/dry spell length diagnostic plot
#'
#' Observed-vs-simulated comparison of average wet and dry spell lengths, faceted by
#' spell type/statistic. Uses dummy points to enforce symmetric axes per facet.
#'
#' @param stats_wetdry Data frame of wet/dry diagnostics including spell statistics.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
#' @import tidyr
.create_spell_length_plot <- function(stats_wetdry, plot_config,
                                     show_title, save_plots, output_dir) {

  data_spells <- stats_wetdry %>%
    dplyr::filter(.data$type == "spells")

  dummy_points <- data_spells %>%
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

  p <- ggplot(data_spells, aes(x = .data$Observed, y = .data$Simulated)) +
    plot_config$theme +
    geom_abline(color = "blue") +
    geom_point(data = dummy_points, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot_config$alpha,
      linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = plot_config$alpha, size = 2) +
    facet_wrap(~ stat, ncol = 2, nrow = 1, scales = "free") +
    labs(x = "Observed", y = "Simulated")

  .export_multipanel_plot(
    p = p,
    filename = "spell_length.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Average dry and wet spell length per month, across all grid cells",
    subtitle = plot_config$subtitle,
    output_dir = output_dir
  )
}


#' Create wet/dry day count diagnostic plot
#'
#' Observed-vs-simulated comparison of average monthly wet and dry day counts, faceted
#' by statistic. Uses dummy points to enforce symmetric axes per facet.
#'
#' @param stats_wetdry Data frame of wet/dry diagnostics including day-count statistics.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
#' @import tidyr
.create_wetdry_days_plot <- function(stats_wetdry, plot_config,
                                    show_title, save_plots, output_dir) {

  data_days <- stats_wetdry %>%
    dplyr::filter(.data$type == "days")

  dummy_points <- data_days %>%
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

  p <- ggplot(data_days, aes(x = .data$Observed, y = .data$Simulated)) +
    plot_config$theme +
    geom_abline(color = "blue") +
    geom_point(data = dummy_points, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot_config$alpha,
      linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = plot_config$alpha, size = 2) +
    facet_wrap(~ stat, ncol = 2, nrow = 1, scales = "free") +
    labs(x = "Observed", y = "Simulated")

  .export_multipanel_plot(
    p = p,
    filename = "wetdry_days_count.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Average number of dry and wet days per month across all grid cells",
    subtitle = plot_config$subtitle,
    output_dir = output_dir
  )
}


#' Create cross-grid correlation diagnostic plot
#'
#' Observed-vs-simulated comparison of cross-grid correlations, faceted by the first
#' variable in each correlation pair. Uses dummy points to enforce symmetric axes.
#'
#' @param stats_crosscor Data frame of cross-grid correlation summaries with columns
#'   \code{Observed}, \code{Simulated}, and \code{variable1}.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
.create_crossgrid_cor_plot <- function(stats_crosscor, plot_config,
                                      show_title, save_plots, output_dir) {

  dummy_points <- generate_symmetric_dummy_points(
    df = stats_crosscor,
    facet_var = "variable1",
    x_col = "Observed",
    y_col = "Simulated"
  )

  p <- ggplot(stats_crosscor, aes(x = .data$Observed, y = .data$Simulated)) +
    plot_config$theme +
    geom_abline(color = "blue") +
    geom_point(data = dummy_points, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot_config$alpha,
      linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = plot_config$alpha, size = 2) +
    facet_wrap(~ variable1, ncol = 2, nrow = 2, scales = "free") +
    labs(x = "Observed", y = "Simulated")

  .export_multipanel_plot(
    p = p,
    filename = "crossgrid_correlations.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Cross-grid correlations",
    subtitle = paste0(plot_config$subtitle, "\nCorrelations calculated over daily series"),
    output_dir = output_dir
  )
}


#' Create inter-variable correlation diagnostic plot
#'
#' Observed-vs-simulated comparison of inter-variable correlations, faceted by variable.
#' Uses dummy points to enforce symmetric axes per facet.
#'
#' @param stats_intercor Data frame of inter-variable correlation summaries with columns
#'   \code{Observed}, \code{Simulated}, and \code{variable}.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
.create_intergrid_cor_plot <- function(stats_intercor, plot_config,
                                      show_title, save_plots, output_dir) {

  dummy_points <- generate_symmetric_dummy_points(
    df = stats_intercor,
    facet_var = "variable",
    x_col = "Observed",
    y_col = "Simulated"
  )

  p <- ggplot(stats_intercor, aes(x = .data$Observed, y = .data$Simulated)) +
    plot_config$theme +
    geom_abline(color = "blue") +
    geom_point(data = dummy_points, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot_config$alpha,
      linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = plot_config$alpha, size = 2) +
    facet_wrap(~ variable, ncol = 3, nrow = 2, scales = "free") +
    labs(x = "Observed", y = "Simulated")

  .export_multipanel_plot(
    p = p,
    filename = "intergrid_correlations.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Inter-variable correlations",
    subtitle = paste0(plot_config$subtitle, "\nCorrelations calculated over daily series"),
    output_dir = output_dir
  )
}


#' Create conditional precipitation correlation diagnostic plot
#'
#' Observed-vs-simulated within-grid correlations between precipitation and other variables,
#' stratified by precipitation regime (e.g., all/wet/dry). Faceted by regime (rows) and
#' variable pair (columns).
#'
#' @param stats_precip_cor_cond Data frame of conditional correlation summaries with columns
#'   \code{variable1}, \code{variable2}, \code{id1}, \code{id2}, \code{regime},
#'   \code{Observed}, and \code{Simulated}.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
.create_precip_cond_cor_plot <- function(stats_precip_cor_cond, plot_config,
                                        show_title, save_plots, output_dir) {

  dat <- stats_precip_cor_cond %>%
    dplyr::mutate(variable = paste0(.data$variable1, ":", .data$variable2)) %>%
    dplyr::filter(.data$id1 == .data$id2)

  dummy_points <- generate_symmetric_dummy_points(
    df = dat,
    facet_var = "variable",
    x_col = "Observed",
    y_col = "Simulated"
  )

  p <- ggplot(dat, aes(x = .data$Observed, y = .data$Simulated)) +
    plot_config$theme +
    geom_abline(color = "blue") +
    geom_point(data = dummy_points, color = "blue", alpha = 0) +
    stat_summary(
      geom = "linerange",
      fun.max = max,
      fun.min = min,
      alpha = plot_config$alpha,
      linewidth = 1.5
    ) +
    stat_summary(fun = "median", geom = "point", alpha = plot_config$alpha, size = 2) +
    facet_grid(regime ~ variable, scales = "free") +
    labs(x = "Observed", y = "Simulated")

  .export_multipanel_plot(
    p = p,
    filename = "precip_conditional_correlations.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Conditional precip-variable correlations (within-grid)",
    subtitle = "Rows: all/wet/dry. Wet uses log1p(precip) if enabled.",
    output_dir = output_dir
  )
}


#' Create monthly pattern plot for a single variable
#'
#' Compares observed and simulated monthly distributions (boxplots) for a given variable,
#' faceted by statistic. Simulated distributions are shown across realizations.
#'
#' @param stats_mon_aavg_sim Data frame of simulated monthly aggregated statistics.
#' @param stats_mon_aavg_obs Data frame of observed monthly aggregated statistics.
#' @param variable Character; variable name to plot (must exist in the inputs).
#' @param plot_config List of plotting configuration options (theme, alpha, colors, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
.create_monthly_pattern_plot <- function(stats_mon_aavg_sim, stats_mon_aavg_obs,
                                        variable, plot_config,
                                        show_title, save_plots, output_dir) {

  dat_sim <- stats_mon_aavg_sim %>%
    dplyr::filter(.data$variable == !!variable) %>%
    dplyr::mutate(
      rlz = as.integer(.data$rlz),
      type = "Simulated"
    ) %>%
    dplyr::rename(value = .data$Simulated)

  dat_obs <- stats_mon_aavg_obs %>%
    dplyr::mutate(rlz = 0L, .before = .data$year) %>%
    dplyr::filter(.data$variable == !!variable, !is.nan(.data$Observed)) %>%
    dplyr::mutate(type = "Observed") %>%
    dplyr::rename(value = .data$Observed)

  dat <- dplyr::bind_rows(dat_sim, dat_obs)

  p <- ggplot(dat, aes(
    x = as.factor(.data$mon),
    y = .data$value,
    fill = .data$type,
    color = .data$type
  )) +
    plot_config$theme +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ stat, scales = "free", ncol = 2, nrow = 2) +
    scale_fill_manual("", values = plot_config$colors) +
    scale_color_manual("", values = plot_config$colors) +
    stat_summary(
      fun = "mean",
      size = 3,
      geom = "point",
      position = position_dodge(0.8),
      shape = 18
    ) +
    labs(x = "", y = "") +
    theme(
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 12)
    ) +
    scale_x_discrete(labels = substr(month.name, 1, 1))

  .export_multipanel_plot(
    p = p,
    filename = paste0("annual_pattern_", variable, ".png"),
    show_title = show_title,
    save_plots = save_plots,
    title = paste0("Monthly patterns for ", variable),
    subtitle = paste0(plot_config$subtitle, "\nResults averaged across all grid cells"),
    output_dir = output_dir
  )
}


#' Create monthly cycle diagnostic plot
#'
#' Plots annual cycles by month for each variable: simulated cycles are shown for all
#' realizations, while the observed cycle is shown as a single reference line.
#'
#' @param daily_stats_season Data frame of daily seasonal statistics including mean values.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title/subtitle.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object (returned invisibly by \code{.export_multipanel_plot()}).
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
.create_monthly_cycle_plot <- function(daily_stats_season, plot_config,
                                      show_title, save_plots, output_dir) {

  sim_avg <- daily_stats_season %>%
    dplyr::group_by(.data$rlz, .data$mon, .data$variable) %>%
    dplyr::filter(.data$stat == "mean") %>%
    dplyr::summarize(value = mean(.data$Simulated, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(type = "Simulated")

  obs_avg <- daily_stats_season %>%
    dplyr::group_by(.data$mon, .data$variable) %>%
    dplyr::filter(.data$stat == "mean") %>%
    dplyr::summarize(value = mean(.data$Observed, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(type = "Observed")

  p <- ggplot(sim_avg, aes(x = as.factor(.data$mon), y = .data$value)) +
    plot_config$theme +
    facet_wrap(~ variable, scales = "free", ncol = 2, nrow = 2) +
    geom_line(aes(group = .data$rlz, color = .data$rlz), alpha = 0.8) +
    geom_line(
      data = obs_avg,
      color = "black",
      group = 1,
      linewidth = 1.25
    ) +
    scale_x_discrete(labels = substr(month.name, 1, 1)) +
    labs(x = "", y = "") +
    guides(color = "none")

  .export_multipanel_plot(
    p = p,
    filename = "monthly_cycle.png",
    show_title = show_title,
    save_plots = save_plots,
    title = "Annual cycles of variables",
    subtitle = paste0(plot_config$subtitle, "\nResults averaged over all grid cells and across each month"),
    output_dir = output_dir
  )
}


#' Create annual mean precipitation diagnostic plot
#'
#' Compares annual mean precipitation time series: simulated realizations are plotted as
#' multiple lines/points and the observed series is overlaid as a single reference series.
#'
#' @param stats_annual_aavg_sim Data frame of simulated annual aggregated statistics.
#' @param stats_annual_aavg_obs Data frame of observed annual aggregated statistics.
#' @param plot_config List of plotting configuration options (theme, alpha, subtitle, etc.).
#' @param show_title Logical; if \code{TRUE}, adds title to the plot.
#' @param save_plots Logical; if \code{TRUE}, writes plot to \code{output_dir}.
#' @param output_dir Character; output directory for saved plots.
#'
#' @return ggplot object returned invisibly.
#'
#' @keywords internal
#' @import ggplot2
#' @import dplyr
.create_annual_precip_plot <- function(stats_annual_aavg_sim, stats_annual_aavg_obs,
                                      plot_config, show_title, save_plots, output_dir) {

  sim_precip <- stats_annual_aavg_sim %>%
    dplyr::filter(.data$stat == "mean", .data$variable == "precip")

  obs_precip <- stats_annual_aavg_obs %>%
    dplyr::filter(.data$stat == "mean", .data$variable == "precip")

  p <- ggplot2::ggplot() +
    plot_config$theme +
    ggplot2::geom_point(
      data = sim_precip,
      ggplot2::aes(x = .data$year, y = .data$Simulated, group = .data$rlz),
      size = 0.5,
      alpha = 0.3
    )

  p <- .add_line_if_possible(
    p = p,
    data = sim_precip,
    mapping = ggplot2::aes(x = .data$year, y = .data$Simulated, group = .data$rlz),
    group_col = "rlz"
  )

  # observed
  p <- p +
    ggplot2::geom_point(
      data = obs_precip,
      ggplot2::aes(x = .data$year, y = .data$Observed),
      size = 0.5,
      alpha = 0.3
    )

  if (nrow(obs_precip) >= 2) {
    p <- p + ggplot2::geom_line(
      data = obs_precip,
      ggplot2::aes(x = .data$year, y = .data$Observed, group = 1),
      linewidth = 1.25
    )
  }

  p <- p + ggplot2::labs(x = "Serial year", y = "mm/day")


  if (show_title) {
    p <- p + labs(title = "Annual mean precipitation")
  }

  if (save_plots && !is.null(output_dir)) {
    ggsave(
      filename = file.path(output_dir, "annual_precip.png"),
      plot = p,
      height = 4,
      width = 8,
      dpi = 300
    )
  }

  invisible(p)
}



#' Generate Dummy Points to Enforce Symmetric Facet Axes
#'
#' @description
#' Creates invisible points with x = y at the min/max range per facet to force ggplot
#' to use symmetric x/y limits within each facet when scales are free.
#'
#' @details
#' The output includes two points per facet (min and max). These points can be added
#' with \code{geom_blank()} to enforce consistent axis limits without affecting
#' the visible data.
#'
#' @param df Data frame containing the plotted data.
#' @param facet_var Character; name of the facet column in \code{df}.
#' @param x_col Character; name of the x column in \code{df}.
#' @param y_col Character; name of the y column in \code{df}.
#'
#' @return Data frame with columns \code{facet_var}, \code{Observed}, \code{Simulated}
#'   (using \code{x_col} and \code{y_col} names) containing min/max dummy points.
#'
#' @examples
#' df <- data.frame(
#'   variable = c("precip", "precip", "temp", "temp"),
#'   Observed = c(1, 5, 10, 12),
#'   Simulated = c(0.5, 6, 9, 13)
#' )
#' generate_symmetric_dummy_points(df, "variable", "Observed", "Simulated")
#'
#' @import dplyr
#' @import rlang
#' @export
generate_symmetric_dummy_points <- function(df, facet_var, x_col, y_col) {

  if (!is.data.frame(df)) stop("df must be a data.frame.", call. = FALSE)
  if (!is.character(facet_var) || length(facet_var) != 1L) stop("facet_var must be a single character string.", call. = FALSE)
  if (!is.character(x_col) || length(x_col) != 1L) stop("x_col must be a single character string.", call. = FALSE)
  if (!is.character(y_col) || length(y_col) != 1L) stop("y_col must be a single character string.", call. = FALSE)

  if (!facet_var %in% names(df)) stop("facet_var not found in df.", call. = FALSE)
  if (!x_col %in% names(df)) stop("x_col not found in df.", call. = FALSE)
  if (!y_col %in% names(df)) stop("y_col not found in df.", call. = FALSE)

  # compute symmetric limits per facet
  out <- df %>%
    dplyr::group_by(.data[[facet_var]]) %>%
    dplyr::summarise(
      .min = min(c(.data[[x_col]], .data[[y_col]]), na.rm = TRUE),
      .max = max(c(.data[[x_col]], .data[[y_col]]), na.rm = TRUE),
      .groups = "drop"
    )

  # two points per facet: (min,min) and (max,max)
  out <- dplyr::bind_rows(
    out %>% dplyr::transmute(!!facet_var := .data[[facet_var]], !!x_col := .data$.min, !!y_col := .data$.min),
    out %>% dplyr::transmute(!!facet_var := .data[[facet_var]], !!x_col := .data$.max, !!y_col := .data$.max)
  )

  out
}

#' Add geom_line only for groups with >= 2 observations
#' @keywords internal
.add_line_if_possible <- function(p, data, mapping, group_col) {

  if (is.null(data) || nrow(data) < 2) return(p)
  if (!group_col %in% names(data)) return(p)

  # count rows per group and keep groups with >=2 points
  tab <- table(data[[group_col]])
  keep <- names(tab)[tab >= 2]

  if (length(keep) == 0) return(p)

  data2 <- data[data[[group_col]] %in% keep, , drop = FALSE]
  if (nrow(data2) < 2) return(p)

  p + ggplot2::geom_line(data = data2, mapping = mapping, inherit.aes = FALSE)
}


