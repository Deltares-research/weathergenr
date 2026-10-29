#' Apply the compact diagnostic presentation profile
#'
#' Restyles returned package plots without regenerating weather or statistics.
#' This explicit opt-in leaves all plotting entry-point defaults unchanged.
#'
#' @param plot A returned weather diagnostic ggplot or patchwork.
#' @param id Stable returned plot ID (not a filename).
#' @param plot_style Partial theme and semantic colors, as documented in
#'   \code{\link{export_weather_plot}}. Explicit theme entries win over this profile.
#' @param panel_layout Optional named arguments for \code{\link{arrange_weather_panels}},
#'   excluding `plots`. NULL uses its defaults for faceted diagnostics. Local
#'   wavelet and WARM summary layouts are preserved.
#' @return A styled ggplot or patchwork with automatic export geometry.
#' @details
#' Observations and agreement guides use `observed` (blue by default).
#' Members and ribbons use `simulated` (gray); central ensemble summaries use
#' black. Significance thresholds and cone-of-influence guides use `reference`
#' (red); significant wavelet regions retain black contours. These semantics
#' apply only to this profile. Boxplots retain medians, quartiles, whiskers and
#' outliers; their diamonds are means. Comparison circles are ensemble medians.
#' No new ensemble summary is computed. Line widths and point sizes are fixed
#' in millimetres, fonts in points. Existing explicit font sizes in `plot_style`
#' override the profile. Use ordinary ggplot additions for further adjustments.
#' Unit labels use Unicode; for vector export choose a Unicode-capable device
#' such as \code{grDevices::cairo_pdf}.
#' @examples
#' # Minimal synthetic annual diagnostic with two supplied series:
#' annual <- data.frame(year = 1:6, Simulated = c(4, 5, 3, 4, 6, 5),
#'   Observed = c(4.2, 5.1, 3.4, 4.3, 5.8, 5.2))
#' p <- ggplot2::ggplot() + theme_weathergenr() +
#'   ggplot2::geom_line(data = annual[c("year", "Simulated")],
#'     ggplot2::aes(year, Simulated)) +
#'   ggplot2::geom_line(data = annual[c("year", "Observed")],
#'     ggplot2::aes(year, Observed))
#' styled <- style_weather_plot(p, "annual_precip")
#' @export
#' @md
style_weather_plot <- function(plot, id, plot_style = NULL, panel_layout = NULL) {
  ids <- c("obs_power_spectra", "warm_annual_precip", "warm_annual_stats",
    "warm_annual_wavelet", "daily_mean", "daily_sd", "spell_length",
    "wetdry_days_count", "crossgrid", "intergrid", "precip_cond_cor",
    "annual_pattern_precip", "annual_pattern_temp", "monthly_cycle",
    "annual_precip")
  if (!inherits(plot, "ggplot") || !is.character(id) || length(id) != 1L ||
      is.na(id) || (!id %in% ids &&
        !grepl("^annual_pattern_.+$", id))) {
    stop("Supply a returned weather diagnostic and its stable plot ID.",
         call. = FALSE)
  }
  if (is.null(plot_style)) {
    plot_style <- list(theme = attr(plot, "weathergenr_style_theme"))
  }
  .validate_plot_style(plot_style)
  if (length(panel_layout) > 0L && (!is.list(panel_layout) ||
      is.null(names(panel_layout)) || anyDuplicated(names(panel_layout)) ||
      any(!names(panel_layout) %in% setdiff(names(formals(arrange_weather_panels)),
                                           "plots")))) {
    stop("panel_layout must contain named arrange_weather_panels arguments.",
         call. = FALSE)
  }
  # Clone layers before editing: ggproto objects have reference semantics.
  clone_plot <- function(p) {
    p$layers <- lapply(p$layers, function(layer) ggplot2::ggproto(NULL, layer))
    p
  }
  original <- if (inherits(plot, "patchwork")) plot else clone_plot(plot)
  style <- .diagnostic_profile()
  base <- list(page = list(raster_dpi = 300, figure_width_mm = 180),
    charts = list(fonts = list(label = 7, tick = 6),
      lines = list(edge = 0.6, line = 1.2),
      rcparams = list("font.family" = "sans", "xtick.major.size" = 2.5)))
  colors <- list(grid = "#e6edef", axis = "#82939b", muted = "#536875",
                 ink = "#17303c")
  role <- function(name, default) .plot_color(plot_style, name, default)
  style$series_roles$observed <- role("observed", "#0D38E0")
  style$series_roles$stochastic_trace <- role("simulated", "#707070")
  style$series_roles$significance <- role("reference", "red")
  style$figures$warm_annual_precip$colors <- list(
    observed = style$series_roles$observed,
    simulated = style$series_roles$stochastic_trace)
  style$figures$obs_power_spectra$colors$significance <- style$series_roles$significance
  style$figures$obs_power_spectra$colors$cone_of_influence <- style$series_roles$significance
  # The wavelet spectrum historically uses blue; retain the profile's blue
  # unless an observed color was explicitly supplied.
  style$figures$obs_power_spectra$colors$observed <- role("observed", "blue")
  style$figures$obs_power_spectra$power_palette <- plot_style$power_palette
  style$panel_layout <- panel_layout
  # Measure with the final text settings rather than adding fonts after layout.
  style$theme <- plot_style$theme
  names_map <- c(crossgrid = "crossgrid_correlations",
    intergrid = "intergrid_correlations",
    precip_cond_cor = "precip_conditional_correlations")
  target <- if (id %in% names(names_map)) names_map[[id]] else id
  result <- suppressMessages(if (id == "obs_power_spectra") {
    .diagnostic_wavelet(original, style$figures$obs_power_spectra, base)
  } else if (id %in% c("annual_precip", "warm_annual_precip")) {
    .diagnostic_annual(original, style$figures$warm_annual_precip,
                   style$annual_time_series, base)
  } else if (id == "warm_annual_stats") {
    .diagnostic_warm_stats(original, style, base, colors)
  } else if (id == "warm_annual_wavelet") {
    .diagnostic_warm_wavelet(original, style, base, colors)
  } else if (startsWith(id, "annual_pattern")) {
    .diagnostic_boxplot(original, target, style, base, colors)
  } else if (id == "monthly_cycle") {
    .diagnostic_monthly(original, style, base, colors)
  } else {
    .diagnostic_comparison(original, target, style, base, colors)
  })
  if (is.null(attr(result, "weathergenr_panel_layout"))) {
    attr(result, "weathergenr_geometry") <- list(family = "wide", dims = NULL,
      size = c(180, 80) / 25.4, header = FALSE)
  }
  # Pin physical type sizes; the caller's partial theme has the last word.
  th <- ggplot2::theme(text = ggplot2::element_text(size = 7, family = "sans"),
    axis.title = ggplot2::element_text(size = 7),
    axis.text = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 6))
  result <- .style_plot(result, list(theme = th))
  result <- .style_plot(result, plot_style)
  layout <- attr(result, "weathergenr_panel_layout")
  if (!is.null(layout) && nrow(layout$cells) == 1L &&
      !isTRUE(panel_layout$single_panel_heading)) {
    result <- .add_plot_theme(result, ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()))
  }
  result
}

.diagnostic_arrange <- function(plots, style, base, result) {
  plots <- lapply(seq_along(plots), function(i) {
    plot <- .add_plot_theme(plots[[i]], style$theme)
    plot
  })
  # Preserve the first-column y title and repeat x titles below each panel.
  # Other columns keep tick values; only a repeated title is removed.
  args <- c(list(plots = plots), style$panel_layout)
  if (is.null(args$repeat_y_title)) args$repeat_y_title <- FALSE
  combined <- do.call(arrange_weather_panels, args)
  attr(combined, "weathergenr_style_theme") <- attr(result, "weathergenr_style_theme")
  attr(combined, "weathergenr_text_size") <- attr(result, "weathergenr_text_size")
  combined
}

.diagnostic_profile <- function() list(series_roles = list(stochastic_trace = "#707070", stochastic_mean = "#000000",
    observed = "#0D38E0", significance = "red"), facet_style = list(
    header_fill = "white", header_text = "#17303c", border = "#8fa5b1",
    border_linewidth_pt = 0.9, spacing_mm = 4L, horizontal_grid_only = TRUE,
    header_alignment = "left", panel_lettering = TRUE), annual_time_series = list(
    horizontal_grid_only = TRUE, colors = list(border = "#aabcc5",
        ticks = "#82939b", text = "#536875", axis_label = "#17303c",
        grid = "#e6edef"), major_tick_step_years = 5L, minor_tick_step_years = 1L,
    minor_tick_length_pt = 1.5, border_linewidth_pt = 0.9), figures = list(
    monthly_cycle = list(simulated_linewidth_mm = 0.25, simulated_alpha = 0.6,
        observed_linewidth_mm = 0.45, x_label = "Month"), monthly_boxplots = list(
        observed = "#0D38E0", simulated_fill = "#707070", simulated_outline = "black",
        fill_alpha = 0.3, linewidth_mm = 0.3, mean_shape = 18L,
        mean_size_mm = 1.25, outlier_size_mm = 1.3, outlier_alpha = 0.25,
        x_label = "Month"), daily_sd = list(x_label = "Observed daily SD",
        y_label = "Simulated daily SD"), daily_mean = list(range_alpha = 0.85,
        range_linewidth_mm = 0.8, median_size_mm = 1.9, median_stroke_mm = 0.3,
        median_fill = "#707070", median_fill_alpha = 0.3, reference_linewidth_mm = 0.3,
        minor_tick_length_pt = 1.5, x_label = "Observed daily mean",
        x_title_per_panel = TRUE, y_label = "Simulated daily mean"),
    obs_power_spectra = list(colors = list(observed = "blue",
        significance = "red", cone_of_influence = "red", significance_contour = "black"),
        spectrum_line_width = "standard", minor_tick_step_years = 1L,
        minor_tick_length_pt = 1.5, align_period_axes = TRUE,
        show_right_period_labels = FALSE, show_right_period_ticks = FALSE,
        panel_gap_mm = 1.5), warm_annual_precip = list(colors = list(
        simulated = "#707070", observed = "#0D38E0"), simulated_linewidth_mm = 0.25,
        simulated_alpha = 0.6, observed_linewidth_mm = 0.45,
        observed_alpha = 1L, x_label = "Serial year", y_label = "mm day\u207B\u00B9"),
    warm_annual_stats = list(x_label = NULL, realization_alpha = 0.5,
        realization_point_size_mm = 1.8, realization_point_outline = "#606060",
        realization_point_stroke_mm = 0.25, realization_dodge_width = 0.3),
    warm_annual_wavelet = list(ensemble_band_alpha = 0.18, x_label = "Period (years)",
        y_label = "mm\u00B2 day\u207B\u00B2", y_limits = list(
            0L, 1.6), y_tick_step = 0.4, minor_tick_step_years = 1L,
        minor_tick_length_pt = 1.5)))


.diagnostic_facet_theme <- function(style, base, colors) {
pt_to_mm <- 25.4 / 72
panels <- style$facet_style
config <- style$figures$daily_mean
facet_theme <- ggplot2::theme(
  text = ggplot2::element_text(size = base$charts$fonts$label,
    family = base$charts$rcparams[["font.family"]]),
  strip.background = ggplot2::element_rect(fill = panels$header_fill, colour = NA),
  strip.text = ggplot2::element_text(colour = panels$header_text,
    face = "plain", size = base$charts$fonts$label, hjust = 0, vjust = 1,
    margin = ggplot2::margin(t = 4, r = 0, b = 4, l = 0)),
  panel.background = ggplot2::element_rect(fill = "white", colour = NA),
  panel.border = ggplot2::element_rect(fill = NA, colour = panels$border,
    linewidth = panels$border_linewidth_pt * pt_to_mm),
  panel.spacing = grid::unit(panels$spacing_mm, "mm"),
  panel.grid.major.x = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  panel.grid.major.y = ggplot2::element_line(colour = colors$grid,
    linewidth = base$charts$lines$edge * pt_to_mm),
  axis.ticks = ggplot2::element_line(colour = colors$axis,
    linewidth = base$charts$lines$edge * pt_to_mm),
  axis.minor.ticks.x.bottom = ggplot2::element_line(colour = colors$axis,
    linewidth = base$charts$lines$edge * pt_to_mm),
  axis.minor.ticks.y.left = ggplot2::element_line(colour = colors$axis,
    linewidth = base$charts$lines$edge * pt_to_mm),
  axis.ticks.length = grid::unit(base$charts$rcparams[["xtick.major.size"]], "pt"),
  axis.minor.ticks.length = grid::unit(config$minor_tick_length_pt, "pt"),
  axis.text = ggplot2::element_text(colour = colors$muted, size = base$charts$fonts$tick),
  axis.title = ggplot2::element_text(colour = colors$ink, size = base$charts$fonts$label)
)

facet_theme
}


.diagnostic_variable_label <- function(value) {
  parts <- strsplit(as.character(value), ":", fixed = TRUE)[[1]]
  names <- c(precip = "Precipitation", temp = "Temperature",
    temp_min = "Minimum temperature", temp_max = "Maximum temperature")
  paste(vapply(parts, function(part) if (part %in% names(names)) names[[part]] else part,
    character(1)), collapse = "\u2013")
}
.diagnostic_axis_labels <- function(target) {
  noun <- switch(target, daily_mean = "daily mean", daily_sd = "daily SD",
    spell_length = "spell length", wetdry_days_count = "day count", "correlation")
  c(x = paste("Observed", noun), y = paste("Simulated", noun))
}
.diagnostic_numeric_equal <- function(expected, measured) {
  stopifnot(nrow(expected) == nrow(measured))
  for (column in intersect(c("x", "y", "ymin", "ymax"), names(expected))) {
    stopifnot(isTRUE(all.equal(expected[[column]], measured[[column]], tolerance = 0)))
  }
}

.diagnostic_comparison <- function(original, target, style, base, colors) {
facet_theme <- .diagnostic_facet_theme(style, base, colors)
config <- style$figures$daily_mean
  stopifnot(inherits(original, "ggplot"), length(original$layers) == 4L)
  old_build <- ggplot2::ggplot_build(original)
  layout <- old_build$layout$layout
  facet_vars <- intersect(c("variable", "variable1", "stat", "regime"), names(layout))
  main_facet <- setdiff(facet_vars, "regime")
  stopifnot(length(main_facet) == 1L)
  reference <- which(vapply(original$layers, function(x) inherits(x$geom, "GeomAbline"), logical(1)))
  ranges <- which(vapply(original$layers, function(x) inherits(x$geom, "GeomLinerange"), logical(1)))
  medians <- which(vapply(original$layers, function(x) inherits(x$stat, "StatSummary") &&
    inherits(x$geom, "GeomPoint"), logical(1)))
  dummy <- setdiff(seq_along(original$layers), c(reference, ranges, medians))
  stopifnot(length(reference) == 1L, length(ranges) == 1L, length(medians) == 1L, length(dummy) == 1L)
  labels <- .diagnostic_axis_labels(target)
  result <- original + facet_theme + ggplot2::labs(x = labels[["x"]], y = labels[["y"]]) +
    ggplot2::scale_x_continuous(guide = ggplot2::guide_axis(minor.ticks = TRUE)) +
    ggplot2::scale_y_continuous(guide = ggplot2::guide_axis(minor.ticks = TRUE))
  result$layers[[ranges]]$aes_params$colour <- style$series_roles$stochastic_trace
  result$layers[[ranges]]$aes_params$alpha <- config$range_alpha
  result$layers[[ranges]]$aes_params$linewidth <- config$range_linewidth_mm
  result$layers[[medians]]$aes_params$shape <- 21
  result$layers[[medians]]$aes_params$colour <- style$series_roles$stochastic_mean
  result$layers[[medians]]$aes_params$fill <- grDevices::adjustcolor(
    config$median_fill, alpha.f = config$median_fill_alpha)
  result$layers[[medians]]$aes_params$stroke <- config$median_stroke_mm
  result$layers[[medians]]$aes_params$size <- config$median_size_mm
  result$layers[[medians]]$aes_params$alpha <- NULL
  result$layers[[reference]]$aes_params$colour <- style$series_roles$observed
  result$layers[[reference]]$aes_params$linewidth <- config$reference_linewidth_mm
  order <- c(dummy, reference, ranges, medians)
  result$layers <- result$layers[order]
  new_build <- ggplot2::ggplot_build(result)
  for (i in seq_along(order)) .diagnostic_numeric_equal(old_build$data[[order[i]]], new_build$data[[i]])
  stopifnot(identical(result$data, original$data))
  panel_plots <- lapply(seq_len(nrow(layout)), function(panel_index) {
    row <- layout[panel_index, , drop = FALSE]
    keep_rows <- function(data) {
      keep <- rep(TRUE, nrow(data))
      for (facet in intersect(facet_vars, names(data))) {
        keep <- keep & as.character(data[[facet]]) == as.character(row[[facet]])
      }
      data[keep, , drop = FALSE]
    }
    value <- as.character(row[[main_facet]])
    title <- .diagnostic_variable_label(value)
    if (target %in% c("daily_mean", "daily_sd")) {
      units <- .diagnostic_units(value)
      if (!is.null(units)) title <- paste0(title, " (", units, ")")
    } else if (target == "spell_length") {
      title <- paste(value, "spells (days)")
    } else if (target == "wetdry_days_count") {
      title <- paste(value, "days (days month\u207b\u00b9)")
    }
    if ("regime" %in% facet_vars) {
      regime <- as.character(row$regime)
      title <- paste(title, paste0(toupper(substr(regime, 1, 1)), substring(regime, 2), " days"), sep = " \u00b7 ")
    }
    title <- paste0(.diagnostic_panel_letter(panel_index), ") ", title)
    panel_plot <- result
    if (nrow(layout) == 1L && target %in% c("daily_mean", "daily_sd")) {
      units <- .diagnostic_units(value)
      if (!is.null(units)) panel_plot <- panel_plot + ggplot2::labs(
        x = paste0(labels[["x"]], " (", units, ")"),
        y = paste0(labels[["y"]], " (", units, ")"))
    }
    panel_plot$data <- keep_rows(result$data)
    panel_plot$layers <- lapply(result$layers, function(layer) {
      copy <- ggplot2::ggproto(NULL, layer)
      if (is.data.frame(copy$data)) copy$data <- keep_rows(copy$data)
      copy
    })
    panel_plot <- panel_plot + ggplot2::facet_wrap(
      stats::reformulate(main_facet), scales = "free", ncol = 1,
      labeller = ggplot2::as_labeller(stats::setNames(title, value)))
    panel_build <- ggplot2::ggplot_build(panel_plot)
    for (i in seq_along(panel_plot$layers)) {
      expected <- new_build$data[[i]]
      expected <- expected[as.integer(expected$PANEL) == as.integer(row$PANEL), , drop = FALSE]
      .diagnostic_numeric_equal(expected, panel_build$data[[i]])
    }
    for (axis in c("x.range", "y.range")) {
      stopifnot(isTRUE(all.equal(new_build$layout$panel_params[[panel_index]][[axis]],
        panel_build$layout$panel_params[[1]][[axis]], tolerance = 0)))
    }
    panel_plot
  })

.diagnostic_arrange(panel_plots, style, base, result)
}


.diagnostic_boxplot <- function(original, target, style, base, colors) {
facet_theme <- .diagnostic_facet_theme(style, base, colors)
equal_summary <- function(a, b) {
  columns <- intersect(c("x", "y", "ymin", "lower", "middle", "upper", "ymax", "outliers", "group"), names(a))
  for (column in columns) stopifnot(isTRUE(all.equal(a[[column]], b[[column]], tolerance = 0)))
}
  stopifnot(inherits(original, "ggplot"))
  old <- ggplot2::ggplot_build(original)
  palette <- c(Observed = style$series_roles$observed, Simulated = style$series_roles$stochastic_trace)
  result <- original + facet_theme + ggplot2::labs(x = "Month", y = NULL) +
    ggplot2::scale_fill_manual(NULL, values = vapply(palette, grDevices::adjustcolor,
      character(1), alpha.f = 0.30)) +
    ggplot2::scale_colour_manual(NULL, values = c(Observed = palette[["Observed"]], Simulated = "black")) +
    ggplot2::scale_x_discrete(labels = substr(month.abb, 1, 1)) +
    ggplot2::theme(legend.position = "none")
  for (layer in result$layers) {
    layer$aes_params$alpha <- NULL
    if (inherits(layer$geom, "GeomBoxplot")) {
      layer$aes_params$linewidth <- 0.30
      if (utils::packageVersion("ggplot2") >= "4.0.0") {
        layer$geom_params$outlier_gp$size <- style$figures$monthly_boxplots$outlier_size_mm
        layer$geom_params$outlier_gp$alpha <- style$figures$monthly_boxplots$outlier_alpha
        layer$geom_params$outlier_gp$stroke <- 0.15
      } else {
        layer$geom_params$outlier.size <- style$figures$monthly_boxplots$outlier_size_mm
        layer$geom_params$outlier.alpha <- style$figures$monthly_boxplots$outlier_alpha
        layer$geom_params$outlier.stroke <- 0.15
      }
    } else if (inherits(layer$geom, "GeomPoint")) {
      layer$aes_params$size <- 1.25
    }
  }
  built <- ggplot2::ggplot_build(result)
  for (i in seq_along(old$data)) equal_summary(old$data[[i]], built$data[[i]])
  layout <- built$layout$layout
  panel_plots <- lapply(seq_len(nrow(layout)), function(i) {
    value <- as.character(layout$stat[i])
    plot <- result
    plot$data <- result$data[as.character(result$data$stat) == value, , drop = FALSE]
    units <- .diagnostic_units(sub("^annual_pattern_", "", target))
    suffix <- if (is.null(units)) "" else paste0(" (", units, ")")
    labels <- c(mean = paste0("Mean", suffix), sd = paste0("Standard deviation", suffix),
      skew = "Skewness", skewness = "Skewness", cor = "Lag-1 correlation", acf = "Lag-1 correlation",
      acf1 = "Lag-1 correlation")
    title <- if (value %in% names(labels)) labels[[value]] else value
    if (nrow(layout) == 1L) plot <- plot + ggplot2::labs(y = title)
    plot <- plot + ggplot2::facet_wrap(~stat, scales = "free", ncol = 1,
      labeller = ggplot2::as_labeller(stats::setNames(paste0(.diagnostic_panel_letter(i), ") ", title), value)))
    single <- ggplot2::ggplot_build(plot)
    for (j in seq_along(built$data)) {
      expected <- built$data[[j]]
      expected <- expected[as.integer(expected$PANEL) == as.integer(layout$PANEL[i]), , drop = FALSE]
      equal_summary(expected, single$data[[j]])
    }
    plot
  })

.diagnostic_arrange(panel_plots, style, base, result)
}


.diagnostic_monthly <- function(original, style, base, colors) {
facet_theme <- .diagnostic_facet_theme(style, base, colors)
stopifnot(inherits(original, "ggplot"), length(original$layers) == 2L)
old <- ggplot2::ggplot_build(original)
config <- style$figures$monthly_cycle
result <- original + facet_theme + ggplot2::labs(x = config$x_label, y = NULL) +
  ggplot2::theme(legend.position = "none")
for (i in seq_along(result$layers)) {
  layer <- result$layers[[i]]
  stopifnot(inherits(layer$geom, "GeomLine"))
  simulated <- "rlz" %in% names(layer$data)
  layer$mapping$colour <- NULL
  layer$aes_params$colour <- if (simulated) style$series_roles$stochastic_trace else style$series_roles$observed
  layer$aes_params$alpha <- if (simulated) config$simulated_alpha else 1
  layer$aes_params$linewidth <- if (simulated) config$simulated_linewidth_mm else config$observed_linewidth_mm
  result$layers[[i]] <- layer
}
equal_lines <- function(a, b) {
  for (column in c("x", "y", "group")) stopifnot(isTRUE(all.equal(a[[column]], b[[column]], tolerance = 0)))
}
built <- ggplot2::ggplot_build(result)
for (i in seq_along(old$data)) equal_lines(old$data[[i]], built$data[[i]])
layout <- built$layout$layout
panel_plots <- lapply(seq_len(nrow(layout)), function(i) {
  value <- as.character(layout$variable[i])
  label <- .diagnostic_variable_label(value)
  units <- .diagnostic_units(value)
  if (!is.null(units)) label <- paste0(label, " (", units, ")")
  keep <- function(data) data[as.character(data$variable) == value, , drop = FALSE]
  plot <- result
  if (nrow(layout) == 1L) plot <- plot + ggplot2::labs(
    y = label)
  plot$data <- keep(result$data)
  plot$layers <- lapply(result$layers, function(layer) {
    copy <- ggplot2::ggproto(NULL, layer)
    copy$data <- keep(layer$data)
    copy
  })
  plot <- plot + ggplot2::facet_wrap(~variable, scales = "free", ncol = 1,
    labeller = ggplot2::as_labeller(stats::setNames(paste0(.diagnostic_panel_letter(i), ") ", label), value)))
  single <- ggplot2::ggplot_build(plot)
  for (j in seq_along(built$data)) {
    expected <- built$data[[j]]
    expected <- expected[as.integer(expected$PANEL) == as.integer(layout$PANEL[i]), , drop = FALSE]
    equal_lines(expected, single$data[[j]])
  }
  for (axis in c("x.range", "y.range")) stopifnot(isTRUE(all.equal(
    built$layout$panel_params[[i]][[axis]], single$layout$panel_params[[1]][[axis]], tolerance = 0)))
  plot
})

.diagnostic_arrange(panel_plots, style, base, result)
}


.diagnostic_annual <- function(original, style, annual_style, base) {
result <- original
for (i in seq_along(result$layers)) {
  data <- result$layers[[i]]$data
  role <- if ("Observed" %in% names(data)) "observed" else
    if (any(c("Simulated", "series", "rlz") %in% names(data))) "simulated" else "observed"
  result$layers[[i]]$aes_params$colour <- style$colors[[role]]
  result$layers[[i]]$aes_params$linewidth <- style[[paste0(role, "_linewidth_mm")]]
  result$layers[[i]]$aes_params$alpha <- style[[paste0(role, "_alpha")]]
}
years <- unlist(lapply(original$layers, function(layer) layer$data$year), use.names = FALSE)
years <- years[is.finite(years)]
if (length(years) > 0L) {
  first <- ceiling(min(years) / annual_style$major_tick_step_years) * annual_style$major_tick_step_years
  major_years <- if (first <= max(years)) seq(first, max(years),
    by = annual_style$major_tick_step_years) else unique(range(years))
  minor_years <- setdiff(seq(min(years), max(years),
    by = annual_style$minor_tick_step_years), major_years)
} else {
  major_years <- ggplot2::waiver()
  minor_years <- ggplot2::waiver()
}
pt_to_mm <- 25.4 / 72
result <- result +
  ggplot2::scale_x_continuous(
    breaks = major_years, minor_breaks = minor_years,
    guide = ggplot2::guide_axis(minor.ticks = TRUE)
  ) +
  ggplot2::labs(x = style$x_label, y = style$y_label) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      colour = annual_style$colors$grid, linewidth = base$charts$lines$edge * pt_to_mm),
    panel.border = ggplot2::element_rect(fill = NA,
      colour = annual_style$colors$border,
      linewidth = annual_style$border_linewidth_pt * pt_to_mm),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(
      colour = annual_style$colors$ticks, linewidth = base$charts$lines$edge * pt_to_mm),
    axis.minor.ticks.x.bottom = ggplot2::element_line(
      colour = annual_style$colors$ticks, linewidth = base$charts$lines$edge * pt_to_mm),
    axis.ticks.length = grid::unit(base$charts$rcparams[["xtick.major.size"]], "pt"),
    axis.minor.ticks.length = grid::unit(annual_style$minor_tick_length_pt, "pt"),
    axis.text = ggplot2::element_text(colour = annual_style$colors$text),
    axis.title = ggplot2::element_text(colour = annual_style$colors$axis_label)
  )

result
}


.diagnostic_warm_stats <- function(original, style, base, colors) {
facet_theme <- .diagnostic_facet_theme(style, base, colors)
panel_labels <- c(
  mean = "a) Mean", sd = "b) Standard\ndeviation",
  tail_low = "c) Dry extremes\n(lower-tail mass)",
  tail_high = "d) Wet extremes\n(upper-tail mass)"
)
result <- original + facet_theme + ggplot2::labs(x = NULL) + ggplot2::facet_wrap(
  ~par, nrow = 1, labeller = ggplot2::as_labeller(panel_labels)
)
stopifnot(length(result$layers) == 2L)
result$layers[[1]]$aes_params$colour <- style$series_roles$observed
result$layers[[1]]$aes_params$linewidth <- 0.30
point_style <- style$figures$warm_annual_stats
result$layers[[2]]$aes_params$shape <- 21
result$layers[[2]]$aes_params$colour <- point_style$realization_point_outline
result$layers[[2]]$aes_params$fill <- grDevices::adjustcolor(
  style$series_roles$stochastic_trace, alpha.f = point_style$realization_alpha
)
result$layers[[2]]$aes_params$alpha <- NULL
result$layers[[2]]$aes_params$size <- point_style$realization_point_size_mm
result$layers[[2]]$aes_params$stroke <- point_style$realization_point_stroke_mm
# Separate realizations consistently across facets without altering their y-values.
result$layers[[2]]$mapping$group <- ggplot2::aes(group = idx)$group
result$layers[[2]]$position <- ggplot2::position_dodge(width = point_style$realization_dodge_width)
stopifnot(identical(result$data, original$data))
point_build <- ggplot2::ggplot_build(result)$data[[2]]
stopifnot(isTRUE(all.equal(sort(point_build$y), sort(result$data$value),
                          check.attributes = FALSE)))
result
}


.diagnostic_warm_wavelet <- function(original, style, base, colors) {
pt_to_mm <- 25.4 / 72
gws <- original
gws$layers[[1]]$aes_params$fill <- style$series_roles$stochastic_trace
gws$layers[[1]]$aes_params$alpha <- style$figures$warm_annual_wavelet$ensemble_band_alpha
gws$layers[[2]]$aes_params$colour <- style$series_roles$stochastic_mean
gws$layers[[3]]$aes_params$colour <- style$series_roles$observed
gws$layers[[4]]$aes_params$colour <- style$series_roles$significance
for (i in 2:4) gws$layers[[i]]$aes_params$linewidth <- c(0.35, 0.45, 0.30)[[i - 1L]]
wavelet_style <- style$figures$warm_annual_wavelet
values <- unlist(lapply(ggplot2::ggplot_build(original)$data, function(data) {
  unlist(data[intersect(c("y", "ymin", "ymax"), names(data))], use.names = FALSE)
}), use.names = FALSE)
maximum <- max(c(0, values[is.finite(values)]))
if (maximum == 0) maximum <- 1
magnitude <- 10^floor(log10(maximum / 4))
steps <- c(1, 2, 4, 5, 10) * magnitude
step <- steps[steps >= maximum / 4][[1]]
wavelet_style$y_tick_step <- step
wavelet_style$y_limits <- c(0, ceiling(maximum / step) * step)
period_range <- range(gws$layers[[2]]$data$period)
major_periods <- ggplot2::ggplot_build(original)$layout$panel_params[[1]]$x$get_breaks()
major_periods <- major_periods[is.finite(major_periods)]
minor_periods <- setdiff(seq(ceiling(period_range[1]), floor(period_range[2]),
  by = wavelet_style$minor_tick_step_years), major_periods)
gws <- gws + ggplot2::scale_x_continuous(
  breaks = major_periods, minor_breaks = minor_periods,
  guide = ggplot2::guide_axis(minor.ticks = TRUE)
) + ggplot2::scale_y_continuous(
    breaks = seq(wavelet_style$y_limits[1], wavelet_style$y_limits[2],
      by = wavelet_style$y_tick_step), expand = c(0, 0)
  ) + ggplot2::coord_cartesian(ylim = wavelet_style$y_limits) +
  ggplot2::labs(x = wavelet_style$x_label, y = wavelet_style$y_label) +
  ggplot2::theme(
    axis.ticks.x = ggplot2::element_line(colour = colors$axis,
      linewidth = base$charts$lines$edge * pt_to_mm),
    axis.minor.ticks.x.bottom = ggplot2::element_line(colour = colors$axis,
      linewidth = base$charts$lines$edge * pt_to_mm),
    axis.ticks.length.x = grid::unit(base$charts$rcparams[["xtick.major.size"]], "pt"),
    axis.minor.ticks.length = grid::unit(wavelet_style$minor_tick_length_pt, "pt")
  )
stopifnot(identical(gws$data, original$data))
stopifnot(isTRUE(all.equal(
  ggplot2::ggplot_build(gws)$layout$panel_params[[1]]$y.range,
  wavelet_style$y_limits, check.attributes = FALSE
)))
gws
}


.diagnostic_wavelet <- function(original, variant, base) {
left <- original[[1]]
if (!is.null(variant$power_palette)) {
  scale <- left$scales$get_scales("fill")
  left <- left + ggplot2::scale_fill_gradientn(colours = variant$power_palette,
    limits = scale$limits, oob = scale$oob, na.value = scale$na.value,
    name = scale$name, breaks = scale$breaks, labels = scale$labels)
}
right <- original[[2]]
left$layers <- lapply(left$layers, function(layer) ggplot2::ggproto(NULL, layer))
right$layers <- lapply(right$layers, function(layer) ggplot2::ggproto(NULL, layer))

left$layers[[2]]$aes_params$colour <- variant$colors$cone_of_influence
left$layers[[3]]$aes_params$colour <- variant$colors$significance_contour
right$layers[[1]]$aes_params$colour <- variant$colors$observed
right$layers[[2]]$aes_params$colour <- variant$colors$significance
stopifnot(variant$spectrum_line_width == "standard")
for (i in 1:2) {
  right$layers[[i]]$aes_params$linewidth <- base$charts$lines$line * 25.4 / 72
}

# The power field includes its outer cell edges. Use those same bounds in the
# spectrum instead of truncating its axis to the wavelet sample centres.
limits <- range(c(left$data$ymin, left$data$ymax))
old_build <- ggplot2::ggplot_build(left)
period_breaks <- as.numeric(old_build$layout$panel_params[[1]]$y$get_labels())
period_breaks <- period_breaks[is.finite(period_breaks)]
time_breaks <- old_build$layout$panel_params[[1]]$x$get_breaks()
time_breaks <- time_breaks[is.finite(time_breaks)]
minor_period <- seq(ceiling(limits[1]), floor(limits[2]), by = variant$minor_tick_step_years)
minor_period <- setdiff(minor_period, period_breaks)
minor_time <- seq(ceiling(min(left$data$xmin)), floor(max(left$data$xmax)),
                  by = variant$minor_tick_step_years)
minor_time <- setdiff(minor_time, time_breaks)
left <- suppressMessages(left +
  ggplot2::scale_y_reverse(breaks = period_breaks, minor_breaks = minor_period,
    expand = c(0, 0), guide = ggplot2::guide_axis(minor.ticks = TRUE)) +
  ggplot2::scale_x_continuous(breaks = time_breaks, minor_breaks = minor_time,
    expand = c(0, 0), guide = ggplot2::guide_axis(minor.ticks = TRUE)))
right <- suppressMessages(right +
  ggplot2::scale_x_reverse(breaks = period_breaks, minor_breaks = minor_period,
    expand = c(0, 0), guide = ggplot2::guide_axis(minor.ticks = TRUE)) +
  ggplot2::coord_flip(xlim = if (isTRUE(variant$align_period_axes)) limits else range(right$data$period),
                      expand = FALSE))
ticks <- ggplot2::theme(
  axis.ticks = ggplot2::element_line(linewidth = 0.25),
  axis.ticks.length = grid::unit(base$charts$rcparams[["xtick.major.size"]], "pt"),
  axis.minor.ticks.x.bottom = ggplot2::element_line(linewidth = 0.25),
  axis.minor.ticks.y.left = ggplot2::element_line(linewidth = 0.25),
  axis.minor.ticks.length = grid::unit(variant$minor_tick_length_pt, "pt")
)
left <- left + ticks
right <- right + ticks
if (!isTRUE(variant$show_right_period_labels)) {
  right <- right + ggplot2::theme(axis.text.y = ggplot2::element_blank())
}
if (!isTRUE(variant$show_right_period_ticks)) {
  right <- right + ggplot2::theme(
    axis.ticks.y = ggplot2::element_blank(),
    axis.minor.ticks.y.left = ggplot2::element_blank(),
    axis.minor.ticks.y.right = ggplot2::element_blank(),
    axis.ticks.length.y = grid::unit(0, "pt"),
    axis.minor.ticks.length.y = grid::unit(0, "pt")
  )
}
half_gap_pt <- variant$panel_gap_mm * 72 / 25.4 / 2
left <- left + ggplot2::theme(
  plot.margin = ggplot2::margin(5.5, half_gap_pt, 5.5, 5.5, unit = "pt")
)
right <- right + ggplot2::theme(
  axis.title.y = ggplot2::element_blank(),
  plot.margin = ggplot2::margin(5.5, 5.5, 5.5, half_gap_pt, unit = "pt")
)
# Use ordinary text for both titles so the superscript does not shift the
# power title's baseline relative to the time title.
right <- right + ggplot2::labs(y = "Power (mm\u00b2)")
horizontal_labels <- ggplot2::theme(
  axis.text.x = ggplot2::element_text(vjust = 1),
  axis.title.x = ggplot2::element_text(vjust = 1)
)
left <- left + horizontal_labels
right <- right + horizontal_labels

# Verify that physical period positions coincide, independently of axis labels.
lb <- ggplot2::ggplot_build(left)
rb <- ggplot2::ggplot_build(right)
for (period in period_breaks) {
  ly <- lb$plot$coordinates$transform(data.frame(x = mean(left$data$x), y = -period),
                                     lb$layout$panel_params[[1]])$y
  ry <- rb$plot$coordinates$transform(data.frame(x = -period, y = mean(right$data$gws)),
                                     rb$layout$panel_params[[1]])$y
  stopifnot(isTRUE(all.equal(ly, ry)))
}
stopifnot(identical(left$data, original[[1]]$data),
          identical(right$data, original[[2]]$data))
result <- left + right + patchwork::plot_layout(widths = c(3, 1.2))
attr(result, "weathergenr_style_theme") <- attr(original, "weathergenr_style_theme")
attr(result, "weathergenr_text_size") <- attr(original, "weathergenr_text_size")
result
}


.diagnostic_units <- function(variable) {
  units <- c(precip = "mm day\u207b\u00b9", temp = "\u00b0C",
    temp_min = "\u00b0C", temp_max = "\u00b0C")
  if (variable %in% names(units)) units[[variable]] else NULL
}

.diagnostic_panel_letter <- function(index) {
  label <- ""
  while (index > 0) {
    index <- index - 1L
    label <- paste0(letters[index %% 26L + 1L], label)
    index <- index %/% 26L
  }
  label
}
