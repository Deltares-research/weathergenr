
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

  plots$precip_cond_cor <- create_precip_cond_cor_plot(
    plot.data$stats.precip_cor_cond,
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
    daily.stats.season =  plot.data$daily.stats.season,
    plot.config = plot.config,
    show.title = show.title,
    save.plots = save.plots,
    output.path = output.path
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

    # force symmetric axes per facet
    ggplot2::geom_point(data = dummy.points, color = "blue", alpha = 0) +

    # overlay summary range + median if you still want it
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

#' Create inter-grid correlation plot conditional on precipitation
#' @keywords internal
create_precip_cond_cor_plot <- function(stats.precip_cor_cond, plot.config,
                                        show.title, save.plots, output.path) {

  dat <- stats.precip_cor_cond %>%
    dplyr::mutate(variable = paste0(.data$variable1, ":", .data$variable2)) %>%
    dplyr::filter(.data$id1 == .data$id2)

  dummy.points <- generate_symmetric_dummy_points(
    df = dat,
    facet.var = "variable",
    x.col = "Observed",
    y.col = "Simulated"
  )

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Observed, y = .data$Simulated)) +
    plot.config$theme +
    ggplot2::geom_abline(color = "blue") +
    ggplot2::geom_point(data = dummy.points, color = "blue", alpha = 0) +
    ggplot2::stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                          alpha = plot.config$alpha, linewidth = 1.5) +
    ggplot2::stat_summary(fun = "median", geom = "point",
                          alpha = plot.config$alpha, size = 2) +
    ggplot2::facet_grid(regime ~ variable, scales = "free") +
    ggplot2::labs(x = "Observed", y = "Simulated")

  export_multipanel_plot(
    p = p,
    filename = "precip_conditional_correlations.png",
    show.title = show.title,
    save.plots = save.plots,
    title = "Conditional precip-variable correlations (within-grid)",
    subtitle = "Rows: all/wet/dry. Wet uses log1p(precip) if enabled.",
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
    subtitle = paste0(plot.config$subtitle, "\nResults averaged over all grid cells and across each month"),
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
