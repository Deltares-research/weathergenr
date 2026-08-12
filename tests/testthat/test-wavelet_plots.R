# Structural assertions on the returned ggplot/patchwork objects, plus a
# `ggplot_build()` call so a broken aesthetic mapping fails here rather than at
# print time. Deliberately not image snapshots: `vdiffr` would add a dependency
# and CI renders on three OSes and three R versions, where SVG output drifts.

# ---- fixture -----------------------------------------------------------------

wavelet_fixture <- function() {
  set.seed(20240101)

  # 64 years with an embedded ~8-year cycle, comfortably above the 16-year
  # minimum enforced by analyze_wavelet_spectrum().
  years <- seq_len(64)
  series <- 100 + 15 * sin(2 * pi * years / 8) + rnorm(64, sd = 4)

  w <- analyze_wavelet_spectrum(series, mode = "complete", seed = 1)
  list(series = series, years = years, w = w)
}

# ---- plot_wavelet_global_spectrum --------------------------------------------

testthat::test_that("plot_wavelet_global_spectrum: builds a ggplot from analysis output", {
  fx <- wavelet_fixture()

  p <- plot_wavelet_global_spectrum(
    period = fx$w$period,
    signif = fx$w$gws_signif,
    obs_power = fx$w$gws
  )

  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_equal(p$labels$x, "Period (years)")
  testthat::expect_equal(nrow(p$data), length(fx$w$period))

  # Observed spectrum + significance threshold.
  testthat::expect_length(p$layers, 2)

  testthat::expect_no_error(ggplot2::ggplot_build(p))
})

testthat::test_that("plot_wavelet_global_spectrum: sim_power adds the ensemble-mean layer", {
  fx <- wavelet_fixture()
  n_period <- length(fx$w$period)

  sim_power <- matrix(
    rep(fx$w$gws, times = 5) * rep(c(0.8, 0.9, 1.0, 1.1, 1.2), each = n_period),
    nrow = n_period
  )

  bare <- plot_wavelet_global_spectrum(
    period = fx$w$period, signif = fx$w$gws_signif, obs_power = fx$w$gws
  )
  withsim <- plot_wavelet_global_spectrum(
    period = fx$w$period, signif = fx$w$gws_signif, obs_power = fx$w$gws,
    sim_power = sim_power
  )

  testthat::expect_length(withsim$layers, length(bare$layers) + 1)
  testthat::expect_no_error(ggplot2::ggplot_build(withsim))

  # The added layer carries the row-wise ensemble mean.
  sim_layer <- withsim$layers[[length(withsim$layers)]]
  testthat::expect_equal(sim_layer$data$sim_mu, rowMeans(sim_power))
})

# ---- plot_wavelet_power ------------------------------------------------------

testthat::test_that("plot_wavelet_power: returns a patchwork of the field and spectrum panels", {
  fx <- wavelet_fixture()

  p <- plot_wavelet_power(
    series = fx$series,
    time = fx$years,
    period = fx$w$period,
    power = fx$w$power,
    gws = fx$w$gws,
    gws_signif = fx$w$gws_signif,
    coi = fx$w$coi,
    signif_mask = fx$w$sigm,
    unit = "mm"
  )

  testthat::expect_s3_class(p, "patchwork")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_wavelet_power: time defaults to the series index", {
  fx <- wavelet_fixture()

  args <- list(
    series = fx$series,
    period = fx$w$period,
    power = fx$w$power,
    gws = fx$w$gws,
    gws_signif = fx$w$gws_signif,
    coi = fx$w$coi,
    signif_mask = fx$w$sigm
  )

  testthat::expect_no_error(do.call(plot_wavelet_power, args))
})

testthat::test_that("plot_wavelet_power: malformed inputs are rejected", {
  fx <- wavelet_fixture()

  base_args <- list(
    series = fx$series,
    time = fx$years,
    period = fx$w$period,
    power = fx$w$power,
    gws = fx$w$gws,
    gws_signif = fx$w$gws_signif,
    coi = fx$w$coi,
    signif_mask = fx$w$sigm
  )

  testthat::expect_error(
    do.call(plot_wavelet_power, utils::modifyList(base_args, list(series = "a"))),
    "'series' must be numeric"
  )

  testthat::expect_error(
    do.call(
      plot_wavelet_power,
      utils::modifyList(base_args, list(power = fx$w$power[, -1, drop = FALSE]))
    ),
    "Invalid 'power' dimensions"
  )

  testthat::expect_error(
    do.call(
      plot_wavelet_power,
      utils::modifyList(base_args, list(signif_mask = as.vector(fx$w$sigm)))
    ),
    "'signif_mask' must be a matrix"
  )
})
