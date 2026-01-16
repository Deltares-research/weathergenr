# tests/testthat/test-plot_filter_diagnostics.R

testthat::skip_if_not_installed("weathergenr")
testthat::skip_if_not_installed("ggplot2")
testthat::skip_if_not_installed("tidyr")
testthat::skip_if_not_installed("dplyr")

.ns <- asNamespace("weathergenr")
testthat::skip_if_not(exists("plot_filter_diagnostics", where = .ns))
testthat::skip_if_not(exists("fill_nearest", where = .ns))

plot_filter_diagnostics <- weathergenr::plot_filter_diagnostics

# ------------------------------------------------------------------------------
# Test data factory
# ------------------------------------------------------------------------------

.make_inputs <- function(n_time = 50L, n_rlz = 20L, n_period = 16L, seed = 1L) {
  set.seed(seed)

  obs_series <- as.numeric(stats::rnorm(n_time, mean = 10, sd = 2))

  sim_series <- matrix(
    stats::rnorm(n_time * n_rlz, mean = 10, sd = 2),
    nrow = n_time,
    ncol = n_rlz
  )

  pool <- seq_len(min(8L, n_rlz))

  rel_diff_mean <- stats::runif(n_rlz, -0.2, 0.2) # relative (fraction)
  rel_diff_sd   <- stats::runif(n_rlz, -0.2, 0.2)

  # Tail metrics need: M_obs_low, M_obs_high scalars; M_sim_low/high vectors (len = n_rlz)
  tail_metrics <- list(
    M_obs_low  = 0.25,
    M_obs_high = 0.30,
    M_sim_low  = stats::runif(n_rlz, 0.15, 0.35),
    M_sim_high = stats::runif(n_rlz, 0.20, 0.40)
  )

  power_period <- sort(stats::runif(n_period, 1, 16))
  power_obs    <- abs(stats::rnorm(n_period, mean = 1.0, sd = 0.2))
  power_signif <- abs(stats::rnorm(n_period, mean = 1.2, sd = 0.2))

  # gws_cache: (n_period x n_rlz)
  gws_cache <- matrix(abs(stats::rnorm(n_period * n_rlz, mean = 1, sd = 0.3)),
                      nrow = n_period, ncol = n_rlz)

  list(
    obs_series = obs_series,
    sim_series = sim_series,
    pool = pool,
    rel_diff_mean = rel_diff_mean,
    rel_diff_sd = rel_diff_sd,
    tail_metrics = tail_metrics,
    power_period = power_period,
    power_obs = power_obs,
    power_signif = power_signif,
    gws_cache = gws_cache
  )
}

# ------------------------------------------------------------------------------
# Output structure tests
# ------------------------------------------------------------------------------

testthat::test_that("returns named list of ggplot objects", {
  x <- .make_inputs()

  res <- plot_filter_diagnostics(
    obs_series     = x$obs_series,
    sim_series     = x$sim_series,
    pool           = x$pool,
    rel_diff_mean  = x$rel_diff_mean,
    rel_diff_sd    = x$rel_diff_sd,
    tail_metrics   = x$tail_metrics,
    power_period   = x$power_period,
    power_obs      = x$power_obs,
    power_signif   = x$power_signif,
    gws_cache      = x$gws_cache,
    wavelet_q      = c(0.05, 0.95)
  )

  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("timeseries", "stats", "wavelet_gws"))
  testthat::expect_s3_class(res$timeseries, "ggplot")
  testthat::expect_s3_class(res$stats, "ggplot")
  testthat::expect_s3_class(res$wavelet_gws, "ggplot")
})

testthat::test_that("plots can be built without error", {
  x <- .make_inputs()

  res <- plot_filter_diagnostics(
    obs_series     = x$obs_series,
    sim_series     = x$sim_series,
    pool           = x$pool,
    rel_diff_mean  = x$rel_diff_mean,
    rel_diff_sd    = x$rel_diff_sd,
    tail_metrics   = x$tail_metrics,
    power_period   = x$power_period,
    power_obs      = x$power_obs,
    power_signif   = x$power_signif,
    gws_cache      = x$gws_cache,
    wavelet_q      = c(0.05, 0.95)
  )

  testthat::expect_no_error(ggplot2::ggplot_build(res$timeseries))
  testthat::expect_no_error(ggplot2::ggplot_build(res$stats))
  testthat::expect_no_error(ggplot2::ggplot_build(res$wavelet_gws))
})

testthat::test_that("timeseries plot contains both obs and sim layers", {
  x <- .make_inputs()

  res <- plot_filter_diagnostics(
    obs_series     = x$obs_series,
    sim_series     = x$sim_series,
    pool           = x$pool,
    rel_diff_mean  = x$rel_diff_mean,
    rel_diff_sd    = x$rel_diff_sd,
    tail_metrics   = x$tail_metrics,
    power_period   = x$power_period,
    power_obs      = x$power_obs,
    power_signif   = x$power_signif,
    gws_cache      = x$gws_cache
  )

  # Expect at least 2 geoms (sim ensemble + obs line)
  testthat::expect_gte(length(res$timeseries$layers), 2L)
})

testthat::test_that("wavelet plot includes ribbon and lines", {
  x <- .make_inputs()

  res <- plot_filter_diagnostics(
    obs_series     = x$obs_series,
    sim_series     = x$sim_series,
    pool           = x$pool,
    rel_diff_mean  = x$rel_diff_mean,
    rel_diff_sd    = x$rel_diff_sd,
    tail_metrics   = x$tail_metrics,
    power_period   = x$power_period,
    power_obs      = x$power_obs,
    power_signif   = x$power_signif,
    gws_cache      = x$gws_cache
  )

  # ribbon + mean line + obs + signif => 4 layers expected (implementation has 4)
  testthat::expect_gte(length(res$wavelet_gws$layers), 4L)
})

# ------------------------------------------------------------------------------
# Input validation tests (errors)
# ------------------------------------------------------------------------------

testthat::test_that("errors if obs_series is not numeric vector", {
  x <- .make_inputs()

  testthat::expect_error(
    plot_filter_diagnostics(
      obs_series     = as.character(x$obs_series),
      sim_series     = x$sim_series,
      pool           = x$pool,
      rel_diff_mean  = x$rel_diff_mean,
      rel_diff_sd    = x$rel_diff_sd,
      tail_metrics   = x$tail_metrics,
      power_period   = x$power_period,
      power_obs      = x$power_obs,
      power_signif   = x$power_signif,
      gws_cache      = x$gws_cache
    ),
    "obs_series must be numeric vector"
  )
})

testthat::test_that("errors if sim_series is not numeric matrix", {
  x <- .make_inputs()

  testthat::expect_error(
    plot_filter_diagnostics(
      obs_series     = x$obs_series,
      sim_series     = as.data.frame(x$sim_series),
      pool           = x$pool,
      rel_diff_mean  = x$rel_diff_mean,
      rel_diff_sd    = x$rel_diff_sd,
      tail_metrics   = x$tail_metrics,
      power_period   = x$power_period,
      power_obs      = x$power_obs,
      power_signif   = x$power_signif,
      gws_cache      = x$gws_cache
    ),
    "sim_series must be numeric matrix"
  )

  testthat::expect_error(
    plot_filter_diagnostics(
      obs_series     = x$obs_series,
      sim_series     = matrix("x", nrow = nrow(x$sim_series), ncol = ncol(x$sim_series)),
      pool           = x$pool,
      rel_diff_mean  = x$rel_diff_mean,
      rel_diff_sd    = x$rel_diff_sd,
      tail_metrics   = x$tail_metrics,
      power_period   = x$power_period,
      power_obs      = x$power_obs,
      power_signif   = x$power_signif,
      gws_cache      = x$gws_cache
    ),
    "sim_series must be numeric matrix"
  )
})

testthat::test_that("errors if sim_series rows do not match length(obs_series)", {
  x <- .make_inputs()

  bad_sim <- x$sim_series[-1, , drop = FALSE]

  testthat::expect_error(
    plot_filter_diagnostics(
      obs_series     = x$obs_series,
      sim_series     = bad_sim,
      pool           = x$pool,
      rel_diff_mean  = x$rel_diff_mean,
      rel_diff_sd    = x$rel_diff_sd,
      tail_metrics   = x$tail_metrics,
      power_period   = x$power_period,
      power_obs      = x$power_obs,
      power_signif   = x$power_signif,
      gws_cache      = x$gws_cache
    ),
    "rows must match length\\(obs_series\\)"
  )
})

testthat::test_that("errors on invalid wavelet_q", {
  x <- .make_inputs()

  testthat::expect_error(
    plot_filter_diagnostics(
      obs_series     = x$obs_series,
      sim_series     = x$sim_series,
      pool           = x$pool,
      rel_diff_mean  = x$rel_diff_mean,
      rel_diff_sd    = x$rel_diff_sd,
      tail_metrics   = x$tail_metrics,
      power_period   = x$power_period,
      power_obs      = x$power_obs,
      power_signif   = x$power_signif,
      gws_cache      = x$gws_cache,
      wavelet_q      = c(0.95, 0.05)
    ),
    "wavelet_q must be two probabilities"
  )

  testthat::expect_error(
    plot_filter_diagnostics(
      obs_series     = x$obs_series,
      sim_series     = x$sim_series,
      pool           = x$pool,
      rel_diff_mean  = x$rel_diff_mean,
      rel_diff_sd    = x$rel_diff_sd,
      tail_metrics   = x$tail_metrics,
      power_period   = x$power_period,
      power_obs      = x$power_obs,
      power_signif   = x$power_signif,
      gws_cache      = x$gws_cache,
      wavelet_q      = c(0, 0.95)
    ),
    "wavelet_q must be two probabilities"
  )

  testthat::expect_error(
    plot_filter_diagnostics(
      obs_series     = x$obs_series,
      sim_series     = x$sim_series,
      pool           = x$pool,
      rel_diff_mean  = x$rel_diff_mean,
      rel_diff_sd    = x$rel_diff_sd,
      tail_metrics   = x$tail_metrics,
      power_period   = x$power_period,
      power_obs      = x$power_obs,
      power_signif   = x$power_signif,
      gws_cache      = x$gws_cache,
      wavelet_q      = c(0.1)
    ),
    "wavelet_q must be two probabilities"
  )
})

# ------------------------------------------------------------------------------
# Behavior tests (basic invariants)
# ------------------------------------------------------------------------------

testthat::test_that("handles tiny observed tail mass without warnings", {
  x <- .make_inputs()

  # Force extreme small observed lower-tail mass; should remain numerically stable
  x$tail_metrics$M_obs_low <- 0

  res <- plot_filter_diagnostics(
    obs_series     = x$obs_series,
    sim_series     = x$sim_series,
    pool           = x$pool,
    rel_diff_mean  = x$rel_diff_mean,
    rel_diff_sd    = x$rel_diff_sd,
    tail_metrics   = x$tail_metrics,
    power_period   = x$power_period,
    power_obs      = x$power_obs,
    power_signif   = x$power_signif,
    gws_cache      = x$gws_cache,
    wavelet_q      = c(0.05, 0.95)
  )

  testthat::expect_s3_class(res$stats, "ggplot")

  testthat::expect_no_warning(ggplot2::ggplot_build(res$stats))
  testthat::expect_no_warning(ggplot2::ggplot_build(res$timeseries))
  testthat::expect_no_warning(ggplot2::ggplot_build(res$wavelet_gws))
})


testthat::test_that("works when sim_series has no column names", {
  x <- .make_inputs()
  colnames(x$sim_series) <- NULL

  res <- plot_filter_diagnostics(
    obs_series     = x$obs_series,
    sim_series     = x$sim_series,
    pool           = x$pool,
    rel_diff_mean  = x$rel_diff_mean,
    rel_diff_sd    = x$rel_diff_sd,
    tail_metrics   = x$tail_metrics,
    power_period   = x$power_period,
    power_obs      = x$power_obs,
    power_signif   = x$power_signif,
    gws_cache      = x$gws_cache
  )

  testthat::expect_s3_class(res$timeseries, "ggplot")
  testthat::expect_no_error(ggplot2::ggplot_build(res$timeseries))
})
