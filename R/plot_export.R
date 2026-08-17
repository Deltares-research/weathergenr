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
#' @details
#' Grouped by geom family, because the value a role needs depends on how it is
#' drawn. Each role carries its whole visual specification -- weight *and* alpha
#' together -- rather than leaving the two to be picked independently at the call
#' site, which is how an ensemble came to compete with the observed series it is
#' meant to sit behind.
#'
#' \code{trace_mono} and \code{trace_hued} are the same conceptual role, an
#' ensemble of realization traces, split because the alpha that reads correctly
#' depends on the trace colour. Monochrome grey traces at 0.5 are nearly as heavy
#' as the observed line drawn over them; hue-scaled pastel traces at 0.2 are
#' close to invisible. One value cannot serve both, and trying made
#' \code{annual_precip} and \code{warm_annual_precip} hard to read.
#'
#' @format Nested named list. Outer names are geom families, inner names are
#'   roles, innermost are the aesthetic values that role sets.
#' @keywords internal
#' @noRd
.PLOT_GEOM <- list(

  line = list(
    # The observed series. Deliberately the heaviest line in any figure: it is
    # the reference every other line is read against.
    observed   = list(linewidth = 1.4, alpha = 1.00),

    # An ensemble of realization traces in a single colour (grey or black).
    trace_mono = list(linewidth = 0.4, alpha = 0.30),

    # An ensemble of realization traces mapped to a colour scale. Needs more
    # alpha than trace_mono to read at all, and can afford it: the hue scale
    # keeps it from reading as one dark mass.
    trace_hued = list(linewidth = 0.5, alpha = 0.70),

    # Fixed guides: 1:1 ablines, zero hlines, significance thresholds, the cone
    # of influence. Present but never competing with the data.
    reference  = list(linewidth = 0.8, alpha = 1.00),

    # A linerange spanning the ensemble at one x position.
    range      = list(linewidth = 1.5, alpha = 0.40)
  ),

  point = list(
    summary = list(size = 2, alpha = 0.40),   # ensemble median marker
    member  = list(size = 3, alpha = 0.30)    # one member: jitter, mean diamond
  ),

  area = list(
    ribbon = list(alpha = 0.20),
    box    = list(alpha = 0.20)
  )
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


#' Reference geometry and type size for base-size scaling
#'
#' The 2 x 2 square figure -- the shape most of the evaluation diagnostics take
#' -- anchors the scale at \code{theme_weathergenr()}'s default 12 pt.
#' @keywords internal
#' @noRd
.BASE_SIZE_REF    <- 12
.BASE_SIZE_REF_IN <- sqrt(8 * 9.1)   # 2 x 2 square, titled
.BASE_SIZE_RANGE  <- c(9, 16)


#' Base font size for a figure of a given geometry
#'
#' @description
#' Scales the type size with the figure's linear extent, so text reads at a
#' consistent size once figures are viewed or placed at a common width.
#'
#' @details
#' A fixed point size makes text physically identical in every figure, which
#' sounds like consistency and is not what a reader sees: these are diagnostics
#' viewed fit-to-window or embedded at a common column width, so a figure gets
#' scaled by its own size before anyone reads it. At a fixed 12 pt the apparent
#' size varied about 2.1-fold across the set -- large on the 8 x 4.5 in figures,
#' small on the 12 x 13.1 in conditional-correlation grid.
#'
#' The measure is \code{sqrt(width * height)} rather than width alone because
#' fitting to a window is constrained by both dimensions; the two figures that
#' looked most different differ far more in height than in width.
#'
#' Clamped because compensation is a correction, not a licence: an unbounded
#' rule would put absurd type on a pathological grid, and the residual mismatch
#' at the two extremes is smaller than the one it removes.
#'
#' @param size Numeric vector \code{c(width, height)} in inches.
#'
#' @return Numeric scalar. Base font size in points.
#'
#' @keywords internal
#' @noRd
.base_size_for <- function(size) {

  extent <- sqrt(max(1e-6, size[[1L]] * size[[2L]]))
  bs <- .BASE_SIZE_REF * (extent / .BASE_SIZE_REF_IN)

  min(max(bs, .BASE_SIZE_RANGE[[1L]]), .BASE_SIZE_RANGE[[2L]])
}


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
#' @param out_dir Character or NULL. Directory to write into. \code{NULL}
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
.export_figure <- function(p, filename, out_dir,
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

  if (!isTRUE(save_plots) || is.null(out_dir)) return(invisible(p))

  # plot_config is the evaluation pipeline's carrier for the public plot_dpi and
  # plot_device arguments, so it wins over the defaults above.
  if (!is.null(plot_config$dpi))    dpi    <- plot_config$dpi
  if (!is.null(plot_config$device)) device <- plot_config$device

  if (is.null(size)) {
    if (is.null(dims)) dims <- .panel_dims(p)
    size <- .figure_size(nrow = dims[["nrow"]], ncol = dims[["ncol"]],
                         family = family, header = header)
  }

  # Type is scaled here rather than in the builders because the geometry it has
  # to track is not known until now. Setting `text` alone is enough: every text
  # element in theme_weathergenr() is rel() and inherits from it, so the whole
  # hierarchy moves together. Adding a partial theme() also merges rather than
  # replaces, which keeps the per-plot overrides -- the inside legend on the
  # monthly-pattern figures, the blanked x axis on the WARM stats strip -- that
  # adding a complete theme here would silently discard.
  p <- p + theme(text = element_text(size = .base_size_for(size)))

  args <- list(
    filename = file.path(out_dir, filename),
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
