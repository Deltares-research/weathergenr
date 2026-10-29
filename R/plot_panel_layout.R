#' Arrange plots at a fixed physical panel size
#'
#' Tries column counts in order, measures the complete canvas, and selects the
#' first grid fitting the width limit. Panels fill rows from left to right.
#' Text and marks are never scaled to fit. Empty trailing cells remain empty.
#'
#' @param plots Nonempty list of single-panel ggplots.
#' @param panel_width_mm,panel_height_mm Positive plotting-area dimensions.
#' @param preferred_columns Ordered distinct integers of at least two.
#' @param max_width_mm Maximum complete canvas width, including decorations.
#' @param dpi Resolution of the temporary measurement device.
#' @param single_panel_heading Whether to retain a single panel's facet heading.
#'   Units should also be present in axis labels when the heading is omitted.
#' @param repeat_y_title Whether to repeat y-axis titles outside the first column.
#' @return A patchwork with millimetre measurements in its
#'   `weathergenr_panel_layout` attribute and stored automatic export geometry.
#' @examples
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
#'   ggplot2::geom_point() + theme_weathergenr(base_size = 7)
#' panels <- arrange_weather_panels(list(p, p))
#' @export
#' @md
arrange_weather_panels <- function(plots, panel_width_mm = 70,
                                   panel_height_mm = 70,
                                   preferred_columns = c(3L, 2L),
                                   max_width_mm = 180, dpi = 300,
                                   single_panel_heading = FALSE,
                                   repeat_y_title = TRUE) {
  positive <- function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x > 0
  if (!is.list(plots) || length(plots) == 0L ||
      !all(vapply(plots, inherits, logical(1), what = "ggplot"))) {
    stop("plots must be a nonempty list of single-panel ggplots.", call. = FALSE)
  }
  for (value in list(panel_width_mm, panel_height_mm, max_width_mm, dpi)) {
    if (!positive(value)) {
      stop("Panel dimensions, maximum width and dpi must be positive.",
           call. = FALSE)
    }
  }
  if (!is.numeric(preferred_columns) || length(preferred_columns) == 0L ||
      anyNA(preferred_columns) || any(!is.finite(preferred_columns)) ||
      any(preferred_columns < 2 | preferred_columns %% 1 != 0) ||
      anyDuplicated(preferred_columns)) {
    stop("preferred_columns must contain distinct integers of at least two.",
         call. = FALSE)
  }
  if (!is.logical(single_panel_heading) || length(single_panel_heading) != 1L ||
      is.na(single_panel_heading)) {
    stop("single_panel_heading must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(repeat_y_title) || length(repeat_y_title) != 1L ||
      is.na(repeat_y_title)) {
    stop("repeat_y_title must be TRUE or FALSE.", call. = FALSE)
  }
  count <- length(plots)
  for (plot in plots) {
    if (inherits(plot, "patchwork") ||
        nrow(ggplot2::ggplot_build(plot)$layout$layout) != 1L) {
      stop("Each supplied plot must contain exactly one panel.", call. = FALSE)
    }
  }
  if (count == 1L && !single_panel_heading) {
    plots[[1]] <- plots[[1]] + ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank())
  }
  # Measure on a private device, restoring the caller's device on exit.
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = max_width_mm, height = 300,
                units = "mm", res = dpi,
                type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  candidates <- if (count == 1L) 1L else
    preferred_columns[preferred_columns <= count]
  measured_columns <- numeric()
  measured_widths <- numeric()
  for (columns in candidates) {
    if (columns * panel_width_mm > max_width_mm + 1e-6) next
    rows <- ceiling(count / columns)
    arranged <- lapply(seq_along(plots), function(i) {
      if (repeat_y_title || (i - 1L) %% columns == 0L) plots[[i]] else
        plots[[i]] + ggplot2::labs(y = NULL)
    })
    combined <- patchwork::wrap_plots(arranged, ncol = columns, byrow = TRUE,
      widths = grid::unit(rep(panel_width_mm, columns), "mm"),
      heights = grid::unit(rep(panel_height_mm, rows), "mm"))
    table <- patchwork::patchworkGrob(combined)
    size <- c(width = grid::convertWidth(sum(table$widths), "mm", TRUE),
              height = grid::convertHeight(sum(table$heights), "mm", TRUE))
    measured_columns <- c(measured_columns, columns)
    measured_widths <- c(measured_widths, size[["width"]])
    if (size[["width"]] > max_width_mm + 1e-6) next
    selected <- .choose_panel_grid(count, measured_columns, measured_widths,
                                   max_width_mm)
    stopifnot(selected[["columns"]] == columns, selected[["rows"]] == rows)
    cells <- table$layout[startsWith(table$layout$name, "panel") &
      table$layout$name != "panel-area", , drop = FALSE]
    if (nrow(cells) != count) {
      stop("Can't identify the supplied plotting areas.", call. = FALSE)
    }
    for (i in seq_len(count)) {
      cell <- cells[i, ]
      actual <- c(grid::convertWidth(sum(table$widths[cell$l:cell$r]),
                                    "mm", TRUE),
                  grid::convertHeight(sum(table$heights[cell$t:cell$b]),
                                     "mm", TRUE))
      if (any(abs(actual - c(panel_width_mm, panel_height_mm)) > 1e-6)) {
        stop("The composition changed the requested panel dimensions.",
             call. = FALSE)
      }
    }
    attr(combined, "weathergenr_panel_layout") <- list(
      columns = columns, rows = rows, panel_width_mm = panel_width_mm,
      panel_height_mm = panel_height_mm, canvas_mm = size,
      cells = data.frame(panel = seq_len(count),
        row = (seq_len(count) - 1L) %/% columns + 1L,
        column = (seq_len(count) - 1L) %% columns + 1L))
    attr(combined, "weathergenr_geometry") <- list(
      family = "square", dims = NULL, size = unname(size) / 25.4,
      header = FALSE)
    return(combined)
  }
  stop("Panels do not fit max_width_mm. Increase the width limit or explicitly reduce panel dimensions; multiple panels never fall back to one column.",
       call. = FALSE)
}

# The planner accepts renderer measurements, so both language adapters can
# exercise identical fixtures independently of their font metrics.
.choose_panel_grid <- function(panel_count, candidate_columns, widths_mm,
                               max_width_mm) {
  if (!is.numeric(panel_count) || length(panel_count) != 1L ||
      !is.finite(panel_count) || panel_count < 1 || panel_count %% 1 != 0 ||
      !is.numeric(widths_mm) || length(widths_mm) != length(candidate_columns) ||
      any(!is.finite(widths_mm)) || any(widths_mm <= 0)) {
    stop("Supply a positive panel count and finite measured candidate widths.",
         call. = FALSE)
  }
  if (!is.numeric(candidate_columns) || any(!is.finite(candidate_columns)) ||
      any(candidate_columns < 1 | candidate_columns %% 1 != 0) ||
      !is.numeric(max_width_mm) || length(max_width_mm) != 1L ||
      !is.finite(max_width_mm) || max_width_mm <= 0) {
    stop("Supply integer columns and a positive width limit.", call. = FALSE)
  }
  eligible <- if (panel_count == 1L) candidate_columns == 1L else
    candidate_columns >= 2L & candidate_columns <= panel_count
  fitting <- which(eligible & widths_mm <= max_width_mm + 1e-6)
  if (length(fitting) == 0L) {
    stop("Panels do not fit max_width_mm.", call. = FALSE)
  }
  columns <- candidate_columns[[fitting[[1L]]]]
  c(columns = columns, rows = ceiling(panel_count / columns))
}
