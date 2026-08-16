# Functions tested (relative paths):
# - R/plot_export.R: theme_weathergenr(), .figure_size(), .panel_dims(),
#   .export_figure()

################################################################################
# .figure_size()
################################################################################

testthat::test_that(".figure_size reproduces the exporter's previous arithmetic", {
  # Before R/plot_export.R existed, .export_multipanel_plot() computed
  #   width <- ncol * 4; height <- nrow * 4 + 0.5
  # The "square" family exists to keep exactly that, so the migration is
  # behaviour-preserving for the seven figures that use it. If these three move,
  # every obs-vs-sim panel silently changed shape.
  testthat::expect_equal(unname(.figure_size(2, 2, "square")), c(8, 8.5))
  testthat::expect_equal(unname(.figure_size(1, 2, "square")), c(8, 4.5))
  testthat::expect_equal(unname(.figure_size(2, 3, "square")), c(12, 8.5))
})

testthat::test_that(".figure_size covers the wide and narrow families", {
  testthat::expect_equal(unname(.figure_size(1, 1, "wide")), c(8, 4.5))

  # warm_annual_stats: four single-column panels in one row. The point of the
  # narrow family is that this totals the 8 in it has always been, where the
  # square family would give 16.
  testthat::expect_equal(unname(.figure_size(1, 4, "narrow")), c(8, 4.5))
})

testthat::test_that(".figure_size adds header room to height only", {
  plain  <- .figure_size(2, 2, "square", header = FALSE)
  titled <- .figure_size(2, 2, "square", header = TRUE)

  testthat::expect_equal(titled[["width"]], plain[["width"]])
  testthat::expect_equal(titled[["height"]] - plain[["height"]], 0.6)
})

testthat::test_that(".figure_size keeps panel size independent of grid shape", {
  # The whole premise of deriving size from the grid: a panel is the same size
  # whether it sits in a 1x1 or a 3x1. A fixed-preset scheme would fail this.
  panel_h <- function(nrow) {
    (.figure_size(nrow, 1, "square")[["height"]] - 0.5) / nrow
  }
  testthat::expect_equal(panel_h(1), panel_h(3))

  panel_w <- function(ncol) .figure_size(1, ncol, "square")[["width"]] / ncol
  testthat::expect_equal(panel_w(1), panel_w(4))
})

testthat::test_that(".figure_size caps a pathological panel count", {
  big <- .figure_size(40, 40, "square", max_in = 30)
  testthat::expect_equal(unname(big), c(30, 30))
})

testthat::test_that(".figure_size rejects an unknown family", {
  testthat::expect_error(.figure_size(1, 1, "tall"), "arg")
})


################################################################################
# .panel_dims()
################################################################################

# One fixture per layout the package actually produces. These are the assertions
# that bind the facet_grid defect: the previous implementation read
# p$facet$params$ncol/nrow, which facet_grid never populates.

make_facet_data <- function() {
  data.frame(
    x = 1:24, y = 1:24,
    g = rep(letters[1:4], 6),
    regime = rep(c("dry", "normal", "wet"), 8),
    variable = rep(c("m", "n"), 12)
  )
}

testthat::test_that(".panel_dims reads an explicit facet_wrap grid", {
  d <- make_facet_data()
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ g, ncol = 2, nrow = 2)

  testthat::expect_equal(.panel_dims(p), c(ncol = 2L, nrow = 2L))
})

testthat::test_that(".panel_dims infers a facet_wrap grid given only nrow", {
  # warm_annual_stats is facet_wrap(~par, nrow = 1) with no ncol, so the number
  # of columns exists only in the built layout.
  d <- make_facet_data()
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ g, nrow = 1)

  testthat::expect_equal(.panel_dims(p), c(ncol = 4L, nrow = 1L))
})

testthat::test_that(".panel_dims reads a facet_grid layout", {
  # The regression test for the reported defect. precip_conditional_correlations
  # is facet_grid(regime ~ variable) with three regimes; the previous exporter
  # returned 2x2 here and wrote three rows of panels into a two-row canvas.
  d <- make_facet_data()
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(regime ~ variable)

  testthat::expect_equal(.panel_dims(p), c(ncol = 2L, nrow = 3L))

  # And the consequence the defect actually had: the figure is taller than the
  # 2x2 fallback would have made it.
  dims <- .panel_dims(p)
  testthat::expect_gt(
    .figure_size(dims[["nrow"]], dims[["ncol"]], "square")[["height"]],
    .figure_size(2, 2, "square")[["height"]]
  )
})

testthat::test_that(".panel_dims returns a single panel for an unfaceted plot", {
  d <- make_facet_data()
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) + ggplot2::geom_point()

  testthat::expect_equal(.panel_dims(p), c(ncol = 1L, nrow = 1L))
})

testthat::test_that(".panel_dims handles a plot built from an empty ggplot()", {
  # .create_annual_precip_plot() starts from a bare ggplot(), whose $data is a
  # waiver rather than a data frame.
  d <- make_facet_data()
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = d, ggplot2::aes(.data$x, .data$y))

  testthat::expect_equal(.panel_dims(p), c(ncol = 1L, nrow = 1L))
})

testthat::test_that(".panel_dims treats a patchwork as a single canvas", {
  # patchwork exposes no stable accessor for its grid, so the contract is that
  # it degrades to 1x1 without erroring and the call site passes `size`.
  d <- make_facet_data()
  one <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) + ggplot2::geom_point()

  pw <- one + one
  testthat::expect_no_error(.panel_dims(pw))
  testthat::expect_equal(.panel_dims(pw), c(ncol = 1L, nrow = 1L))
})


################################################################################
# theme_weathergenr()
################################################################################

testthat::test_that("theme_weathergenr is a complete theme distinguishable from its neighbours", {
  th <- theme_weathergenr()

  testthat::expect_s3_class(th, "theme")

  # panel.border$colour is the discriminator the plot tests use, because it
  # separates all three themes this package used to mix. panel.background$fill
  # does not: theme_light() and theme_bw() both give white.
  testthat::expect_identical(th$panel.border$colour,
                             ggplot2::theme_light()$panel.border$colour)
  testthat::expect_false(
    identical(th$panel.border$colour, ggplot2::theme_bw()$panel.border$colour)
  )
  testthat::expect_false(
    identical(th$panel.border$colour, ggplot2::theme_grey()$panel.border$colour)
  )
})

testthat::test_that("theme_weathergenr defaults to base_size 12 and scales with it", {
  # 12 matches what the evaluation pipeline produced before this existed, so the
  # standardization moves geometry without also moving type size.
  testthat::expect_equal(theme_weathergenr()$text$size, 12)
  testthat::expect_equal(theme_weathergenr(base_size = 9)$text$size, 9)
})

testthat::test_that("theme_weathergenr title sizes are relative, not absolute", {
  # They were flat 14 and 10 points, which would have been stranded by any
  # change to base_size.
  th <- theme_weathergenr()
  testthat::expect_s3_class(th$plot.title$size, "rel")
  testthat::expect_s3_class(th$plot.subtitle$size, "rel")
})

testthat::test_that("theme_weathergenr rejects a nonsense base_size", {
  testthat::expect_error(theme_weathergenr(base_size = -1), "positive")
  testthat::expect_error(theme_weathergenr(base_size = c(10, 12)), "single")
})


################################################################################
# .export_figure()
################################################################################

# A recorder, not a counter. The previous export test counted calls, and because
# the helper called ggplot2::ggsave -- namespace-qualified, so invisible to
# local_mocked_bindings -- the count it reported came from a single call site.

local_ggsave_recorder <- function(env = parent.frame()) {
  calls <- new.env(parent = emptyenv())
  calls$args <- list()

  testthat::local_mocked_bindings(
    ggsave = function(...) {
      a <- list(...)
      calls$args[[length(calls$args) + 1L]] <- a
      invisible(NULL)
    },
    .package = "weathergenr",
    .env = env
  )

  calls
}

testthat::test_that(".export_figure always sets units, dpi and background", {
  rec <- local_ggsave_recorder()
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()

  .export_figure(p, "f.png", output_dir = tempdir(), family = "wide")

  testthat::expect_length(rec$args, 1L)
  a <- rec$args[[1]]
  testthat::expect_identical(a$units, "in")
  testthat::expect_identical(a$bg, "white")
  testthat::expect_equal(a$dpi, 300)
  testthat::expect_equal(a$width, 8)
  testthat::expect_equal(a$height, 4.5)
})

testthat::test_that(".export_figure omits device rather than passing NULL", {
  # ggsave() infers the device from the extension when the argument is absent;
  # an explicit NULL is not the same thing.
  rec <- local_ggsave_recorder()
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()

  .export_figure(p, "f.png", output_dir = tempdir(), device = NULL)
  testthat::expect_false("device" %in% names(rec$args[[1]]))

  .export_figure(p, "g.png", output_dir = tempdir(), device = grDevices::png)
  testthat::expect_true("device" %in% names(rec$args[[2]]))
})

testthat::test_that(".export_figure lets plot_config carry dpi and device", {
  # This is how the public plot_dpi / plot_device arguments reach the writer.
  rec <- local_ggsave_recorder()
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()

  .export_figure(p, "f.png", output_dir = tempdir(),
                 dpi = 300, plot_config = list(dpi = 150))

  testthat::expect_equal(rec$args[[1]]$dpi, 150)
})

testthat::test_that(".export_figure honours the dims and size overrides", {
  rec <- local_ggsave_recorder()
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()

  # dims overrides detection; the plot itself is unfaceted.
  .export_figure(p, "a.png", output_dir = tempdir(), family = "square",
                 dims = c(ncol = 3L, nrow = 2L))
  testthat::expect_equal(rec$args[[1]]$width, 12)
  testthat::expect_equal(rec$args[[1]]$height, 8.5)

  # size overrides everything, which is what the patchwork call site uses.
  .export_figure(p, "b.png", output_dir = tempdir(), family = "square",
                 dims = c(ncol = 3L, nrow = 2L), size = c(9, 4.5))
  testthat::expect_equal(rec$args[[2]]$width, 9)
  testthat::expect_equal(rec$args[[2]]$height, 4.5)
})

testthat::test_that(".export_figure adds the title block and its height together", {
  rec <- local_ggsave_recorder()
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()

  out <- .export_figure(p, "f.png", output_dir = tempdir(), family = "wide",
                        show_title = TRUE, title = "T", subtitle = "S")

  testthat::expect_identical(out$labels$title, "T")
  testthat::expect_identical(out$labels$subtitle, "S")
  testthat::expect_equal(rec$args[[1]]$height, 5.1)

  # show_title = TRUE with no title is not a title: neither label nor room.
  .export_figure(p, "g.png", output_dir = tempdir(), family = "wide",
                 show_title = TRUE, title = NULL)
  testthat::expect_equal(rec$args[[2]]$height, 4.5)
})

testthat::test_that(".export_figure writes nothing when told not to", {
  rec <- local_ggsave_recorder()
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()

  out1 <- .export_figure(p, "f.png", output_dir = tempdir(), save_plots = FALSE)
  out2 <- .export_figure(p, "f.png", output_dir = NULL)

  testthat::expect_length(rec$args, 0L)
  testthat::expect_s3_class(out1, "ggplot")
  testthat::expect_s3_class(out2, "ggplot")
})
