# ==============================================================================
# Shared figure style and export
# ==============================================================================
#
# The one place ggsave() is called and the one place the house theme is built.
#
# This is not a `*_plots.R` sibling because it belongs to no module: three plot
# modules (evaluate_generator_plots.R, warm_filtering_plots.R, wavelet_plots.R)
# plus the writer for generator.R consume it, and putting it in any one sibling
# would make the other two depend on an unrelated module's file. Centralizing it
# is also what lets generator.R hold no ggsave() and evaluate_generator.R hold no
# ggplot2 object -- both are computational files that previously did.

# ------------------------------------------------------------------------------
# Theme
# ------------------------------------------------------------------------------

#' weathergenr house theme
#'
#' @description
#' A \code{\link[ggplot2]{theme_light}} wrapper applied to every figure this
#' package draws, so the PNGs a run writes share one look.
#'
#' \code{evaluate_weather_generator()} returns its diagnostic plots as ggplot
#' objects, so this is exported: without it, a caller re-rendering one of those
#' plots cannot reproduce the package's styling.
#'
#' @details
#' \code{base_size} is 12 rather than \code{theme_light()}'s native 11 because 12
#' is what the evaluation pipeline has always produced. The title and subtitle
#' sizes are \code{\link[ggplot2]{rel}} rather than absolute points so they track
#' \code{base_size}; they were previously fixed at 14 and 10 and would have been
#' stranded by any change to it.
#'
#' Legend text is deliberately not set here. Exactly one exported figure has a
#' visible legend, so its styling stays local to that plot rather than becoming a
#' package-wide rule inferred from a single case.
#'
#' @param base_size Numeric. Base font size in points. Default 12.
#' @param base_family Character. Base font family. Default \code{""}, the
#'   device default.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_weathergenr()
#'
#' @export
#' @import ggplot2
theme_weathergenr <- function(base_size = 12, base_family = "") {

  if (!is.numeric(base_size) || length(base_size) != 1L || !is.finite(base_size) ||
      base_size <= 0) {
    stop("base_size must be a single positive number.", call. = FALSE)
  }

  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(size = rel(1.2)),
      plot.subtitle = element_text(size = rel(0.85))
    )
}


# ------------------------------------------------------------------------------
# Shared geom constants
# ------------------------------------------------------------------------------

#' Shared geom weights, sizes and alphas
#'
#' @description
#' One entry per visual *role*, not per call site, so that -- for example -- the
#' observed reference series is drawn at the same weight in every figure. Before
#' this existed the observed series was 0.9 in one module, 1.25 in another, and
#' the device default in a third.
#'
#' @format Named list of numeric scalars.
#' @keywords internal
#' @noRd
.PLOT_GEOM <- list(
  lwd_observed  = 1.0,   # the observed / reference series
  lwd_simulated = 0.5,   # a single simulated member
  lwd_summary   = 1.5,   # stat_summary linerange spanning the ensemble
  lwd_reference = 0.8,   # abline, hline, significance threshold
  pt_summary    = 2,     # ensemble median marker
  pt_member     = 3,     # single-member marker (jitter point, mean diamond)
  alpha_faint   = 0.2,   # area fills: ribbons, boxplot interiors
  # A bundle of realization lines drawn behind a reference series. Separate from
  # alpha_faint because an area fill and a line bundle do not read the same at
  # the same alpha: 0.2 is right for a ribbon but leaves hue-scaled lines almost
  # invisible, which is what happened to monthly_cycle. Held below the 0.8 the
  # bundles used to carry, where 20 realizations saturate to solid.
  alpha_bundle  = 0.5,
  alpha_member  = 0.3,   # individual member points
  alpha_summary = 0.4    # ensemble range and median summaries
)


# ------------------------------------------------------------------------------
# Figure geometry
# ------------------------------------------------------------------------------

#' Panel geometry by figure family
#'
#' @description
#' Panel size is a property of what is *inside* a panel, so it does not vary with
#' the shape of the grid: an n x m grid is \code{m * panel_w} by
#' \code{n * panel_h} plus a fixed margin allowance for axes and strips.
#'
#' \code{square} reproduces the arithmetic the evaluation exporter used before
#' this file existed (\code{width <- ncol * 4}, \code{height <- nrow * 4 + 0.5}).
#' That is deliberate: those panels plot Observed against Simulated with a 1:1
#' \code{geom_abline()}, which is only readable at an aspect ratio of 1, and
#' keeping the number means any size change that does appear is one somebody
#' chose rather than a side effect.
#'
#' \code{narrow} exists for one figure -- \code{warm_annual_stats}, four panels
#' each holding a single jittered column with the x axis blanked. Under
#' \code{square} it would render 16 inches wide for four columns of dots.
#'
#' @format Named list of lists, each with \code{panel_w}, \code{panel_h},
#'   \code{margin_w} and \code{margin_h} in inches.
#' @keywords internal
#' @noRd
.PANEL_SIZES <- list(
  square = list(panel_w = 4.0, panel_h = 4.0, margin_w = 0.0, margin_h = 0.5),
  wide   = list(panel_w = 8.0, panel_h = 4.0, margin_w = 0.0, margin_h = 0.5),
  narrow = list(panel_w = 2.0, panel_h = 4.0, margin_w = 0.0, margin_h = 0.5)
)

#' Extra height, in inches, for a title and subtitle block
#'
#' The margin above is unconditional and was sized for neither. `show_title` is
#' the pipeline default and three of the subtitles wrap to two lines, so a titled
#' figure needs the room; an untitled one does not.
#' @keywords internal
#' @noRd
.HEADER_ALLOWANCE_IN <- 0.6


#' Figure dimensions for a panel grid
#'
#' @param nrow,ncol Integer. Panel grid dimensions.
#' @param family Character. One of \code{"square"}, \code{"wide"},
#'   \code{"narrow"}; see \code{.PANEL_SIZES}.
#' @param header Logical. If \code{TRUE}, add room for a title/subtitle block.
#' @param max_in Numeric. Upper bound on either dimension, in inches.
#'
#' @return Named numeric vector \code{c(width, height)} in inches.
#'
#' @keywords internal
#' @noRd
.figure_size <- function(nrow = 1L, ncol = 1L,
                         family = c("square", "wide", "narrow"),
                         header = FALSE,
                         max_in = 30) {

  family <- match.arg(family)
  g <- .PANEL_SIZES[[family]]

  nrow <- max(1L, as.integer(nrow))
  ncol <- max(1L, as.integer(ncol))

  width  <- ncol * g$panel_w + g$margin_w
  height <- nrow * g$panel_h + g$margin_h +
    if (isTRUE(header)) .HEADER_ALLOWANCE_IN else 0

  # ggsave() errors above 50 in. Cap well below it so a pathological panel count
  # degrades to a squashed figure rather than to a failed run.
  c(width = min(width, max_in), height = min(height, max_in))
}


#' Infer a plot's panel grid
#'
#' @description
#' Returns the number of panel columns and rows a plot will draw.
#'
#' @details
#' The grid comes from the *built* layout rather than from the facet
#' specification. Reading \code{p$facet$params$ncol} / \code{$nrow} instead is
#' what made \code{facet_grid} figures fall through to a hardcoded 2x2: those
#' fields are populated by \code{facet_wrap} only, so a three-row
#' \code{facet_grid} was being written into a two-row canvas. One
#' \code{ggplot_build()} call resolves \code{facet_wrap} (with or without
#' explicit \code{ncol}/\code{nrow}), \code{facet_grid}, and the unfaceted case
#' alike, and costs nothing next to the render that follows it.
#'
#' @param p A ggplot or patchwork object.
#'
#' @return Named integer vector \code{c(ncol, nrow)}.
#'
#' @keywords internal
#' @noRd
.panel_dims <- function(p) {

  # A patchwork governs its own internal split through plot_layout(), so the
  # outer canvas does not need to know how many plots are inside it. This needs
  # its own branch because a patchwork is also class "ggplot", and because
  # patchwork exposes no stable accessor for its grid: attr(p, "patches") is
  # NULL, p@patches errors, and p[["patches"]] raises "Patchworks can only be
  # indexed with numeric indices". Those call sites pass `size` explicitly.
  if (inherits(p, "patchwork")) return(c(ncol = 1L, nrow = 1L))

  tryCatch({
    lay <- ggplot_build(p)$layout$layout
    c(ncol = max(as.integer(lay$COL)), nrow = max(as.integer(lay$ROW)))
  }, error = function(e) c(ncol = 1L, nrow = 1L))
}


# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

#' Write a figure at the package's standard geometry
#'
#' @description
#' The single \code{ggsave()} call site in the package. Attaches the title block
#' when one is wanted, derives the figure size from the panel grid unless told
#' otherwise, and writes with \code{units}, \code{dpi} and \code{bg} always set
#' explicitly rather than left to \code{ggsave()}'s defaults.
#'
#' @param p A ggplot or patchwork object.
#' @param filename Character. Base filename, including extension.
#' @param output_dir Character or NULL. Directory to write into. \code{NULL}
#'   writes nothing.
#' @param save_plots Logical. \code{FALSE} writes nothing.
#' @param show_title Logical. If \code{TRUE} and \code{title} is not \code{NULL},
#'   attach the title block and allow room for it.
#' @param title,subtitle Character or NULL. Title block contents.
#' @param family Character. Figure family; see \code{.PANEL_SIZES}.
#' @param dims Named integer vector \code{c(ncol, nrow)} or NULL. Overrides
#'   grid detection.
#' @param size Numeric vector \code{c(width, height)} in inches, or NULL.
#'   Overrides both detection and \code{family}.
#' @param dpi Numeric. Raster resolution.
#' @param device Function or NULL. Graphics device passed to
#'   \code{ggplot2::ggsave()}.
#' @param plot_config List or NULL. When it carries \code{dpi} or \code{device},
#'   those win over the arguments above -- this is how the public
#'   \code{plot_dpi} / \code{plot_device} settings reach the writer.
#'
#' @return \code{p}, invisibly.
#'
#' @keywords internal
#' @noRd
#' @import ggplot2
.export_figure <- function(p, filename, output_dir,
                           save_plots = TRUE,
                           show_title = FALSE, title = NULL, subtitle = NULL,
                           family = c("square", "wide", "narrow"),
                           dims = NULL,
                           size = NULL,
                           dpi = 300, device = NULL,
                           plot_config = NULL) {

  family <- match.arg(family)

  header <- isTRUE(show_title) && !is.null(title)
  if (header) p <- p + labs(title = title, subtitle = subtitle)

  if (!isTRUE(save_plots) || is.null(output_dir)) return(invisible(p))

  # plot_config is the evaluation pipeline's carrier for the public plot_dpi and
  # plot_device arguments, so it wins over the defaults above.
  if (!is.null(plot_config$dpi))    dpi    <- plot_config$dpi
  if (!is.null(plot_config$device)) device <- plot_config$device

  if (is.null(size)) {
    if (is.null(dims)) dims <- .panel_dims(p)
    size <- .figure_size(nrow = dims[["nrow"]], ncol = dims[["ncol"]],
                         family = family, header = header)
  }

  args <- list(
    filename = file.path(output_dir, filename),
    plot     = p,
    width    = unname(size[[1L]]),
    height   = unname(size[[2L]]),
    units    = "in",
    dpi      = dpi,
    bg       = "white"
  )

  # Only pass `device` when set: ggsave() infers it from the extension
  # otherwise, and an explicit NULL is not the same as omitting it.
  if (!is.null(device)) args$device <- device

  # Bare `ggsave` (via import(ggplot2)), never ggplot2::ggsave. The qualified
  # form bypasses testthat's local_mocked_bindings(.package = "weathergenr"),
  # which is how the export test came to be observing one of six call sites
  # while the other five went unchecked.
  do.call(ggsave, args)

  invisible(p)
}
