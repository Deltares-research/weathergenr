# Functions tested (relative paths):
# - R/wavelet.R: analyze_wavelet_spectrum(), fill_nearest(), extract_signif_curve(),
#   gws_regrid(), morlet_wavelet(), morlet_parameters(), extract_wavelet_components(),
#   simulate_warm()
# - R/wavelet_plots.R: plot_wavelet_power(), plot_wavelet_global_spectrum()

testthat::test_that("analyze_wavelet_spectrum returns expected structure", {
  series <- sin(seq(0, 4 * pi, length.out = 64))

  out_fast <- analyze_wavelet_spectrum(series, noise = "white", mode = "fast")
  testthat::expect_true(all(c("gws", "period", "power", "coi") %in% names(out_fast)))
  testthat::expect_equal(nrow(out_fast$power), length(out_fast$period))
  testthat::expect_equal(ncol(out_fast$power), length(series))
  testthat::expect_equal(length(out_fast$coi), length(series))

  out_complete <- analyze_wavelet_spectrum(series, noise = "white", mode = "complete", diagnostics = TRUE)
  testthat::expect_true(all(c("power_signif_coi", "wave", "comps", "neff") %in% names(out_complete)))
  testthat::expect_true(is.matrix(out_complete$comps))
  testthat::expect_true(is.list(out_complete$diagnostics))
})

testthat::test_that("fill_nearest fills missing values deterministically", {
  x <- c(NA, 1, NA, 2, NA)
  out <- fill_nearest(x)
  testthat::expect_identical(out, c(1, 1, 1, 2, 2))

  all_na <- c(NA_real_, NA_real_)
  testthat::expect_identical(fill_nearest(all_na), all_na)
})

testthat::test_that("extract_signif_curve uses known keys and handles missing", {
  wavelet_a <- list(gws_signif = c(1, 2, 3))
  testthat::expect_identical(extract_signif_curve(wavelet_a), c(1, 2, 3))

  wavelet_b <- list(signif = c(0.1, 0.2))
  testthat::expect_identical(extract_signif_curve(wavelet_b), c(0.1, 0.2))

  testthat::expect_true(is.null(extract_signif_curve(list())))
})

testthat::test_that("gws_regrid interpolates and respects unmasked option", {
  wavelet <- list(
    period = c(1, 2, 3),
    gws = c(10, 20, 30),
    gws_unmasked = c(11, 21, 31)
  )

  target <- c(1, 1.5, 3)
  out <- gws_regrid(wavelet, target)
  expected <- stats::approx(x = wavelet$period, y = wavelet$gws, xout = target, rule = 2)$y
  testthat::expect_equal(out, expected)

  out_unmasked <- gws_regrid(wavelet, target, use_unmasked = TRUE)
  expected_unmasked <- stats::approx(x = wavelet$period, y = wavelet$gws_unmasked,
                                     xout = target, rule = 2)$y
  testthat::expect_equal(out_unmasked, expected_unmasked)
})

testthat::test_that("morlet_wavelet and morlet_parameters return valid outputs", {
  k <- c(0:4, -3:-1)
  out <- morlet_wavelet(k, scale = 2, k0 = 6)
  testthat::expect_equal(length(out), length(k))
  testthat::expect_true(all(out[k <= 0] == 0))
  testthat::expect_true(any(out[k > 0] != 0))

  params <- morlet_parameters(k0 = 6)
  testthat::expect_true(all(c("fourier_factor", "coi", "dofmin") %in% names(params)))
  testthat::expect_equal(unname(params["dofmin"]), 2)
})

testthat::test_that("extract_wavelet_components reconstructs components and attributes", {
  wave <- matrix(c(1, 2, 3, 4,
                   2, 3, 4, 5,
                   1, 1, 1, 1),
                 nrow = 3, byrow = TRUE)
  scale <- c(1, 2, 3)

  comps <- extract_wavelet_components(
    wave = wave,
    signif_periods = c(1, 3),
    scale = scale,
    series_sd = 2,
    series_mean = 0
  )

  testthat::expect_equal(nrow(comps), ncol(wave))
  testthat::expect_equal(ncol(comps), 3)
  testthat::expect_true(all(c("Component_1", "Component_2", "Noise") %in% colnames(comps)))
  testthat::expect_equal(attr(comps, "n_significant"), 2)

  testthat::expect_warning(
    extract_wavelet_components(
      wave = wave,
      signif_periods = c(1, 1, 3),
      scale = scale,
      series_sd = 2,
      series_mean = 0
    ),
    "duplicates"
  )
})

################################################################################


testthat::test_that("simulate_warm returns deterministic simulations", {
  components <- cbind(
    rep(c(0, 1), length.out = 20),
    rep(c(1, 0), length.out = 20)
  )
  out1 <- simulate_warm(components = components, n = 20, n_sim = 3, seed = 10, verbose = FALSE)
  out2 <- simulate_warm(components = components, n = 20, n_sim = 3, seed = 10, verbose = FALSE)

  testthat::expect_equal(dim(out1), c(20, 3))
  testthat::expect_identical(out1, out2)
})

test_that("simulate_warm bypass mode uses block bootstrap when ARIMA not viable", {
  series_obs <- rnorm(20)
  n <- length(series_obs)

  calls_viable <- 0L
  calls_boot <- 0L

  testthat::local_mocked_bindings(
    .fit_warm_arima_forecast = function(x, max_p, max_q) {
      calls_viable <<- calls_viable + 1L
      NULL
    },
    .block_bootstrap = function(x, n, block_len) {
      calls_boot <<- calls_boot + 1L
      rep(mean(x), n)
    },
    .default_block_len = function(n) 5L,
    .env = asNamespace("weathergenr")
  )

  out <- weathergenr::simulate_warm(
    components = NULL,
    n = n,
    n_sim = 5,
    seed = 1,
    series_obs = series_obs,
    bypass_n = 25L,          # force bypass
    verbose = FALSE
  )

  expect_equal(calls_viable, 1L)
  expect_equal(calls_boot, 5L)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(n, 5L))
})


testthat::test_that("simulate_warm component mode uses block bootstrap when ARIMA not viable, without warnings", {

  skip_if_not_installed("forecast")

  n <- 60L
  # low uniqueness -> not viable in our gating logic
  d1 <- rep(c(0, 1), length.out = n)
  comps <- cbind(D1 = d1)

  calls_viable <- 0L
  calls_fit <- 0L
  calls_boot <- 0L

  testthat::local_mocked_bindings(
    .warm_arima_viable = function(x, min_n, sd_eps, min_unique) {
      calls_viable <<- calls_viable + 1L
      FALSE
    },
    .warm_fit_arima_safe = function(...) {
      calls_fit <<- calls_fit + 1L
      stop("should not be called")
    },
    .warm_block_bootstrap = function(x, n) {
      calls_boot <<- calls_boot + 1L
      # deterministic: return centered pattern (preserves dependence in real impl)
      x[seq_len(n)]
    },
    .env = asNamespace("weathergenr")
  )

  out <- testthat::expect_silent(
    simulate_warm(
      components = comps,
      n = n,
      n_sim = 4L,
      seed = 1L,
      bypass_n = 25L,
      verbose = FALSE
    )
  )

  testthat::expect_equal(dim(out), c(n, 4L))
  testthat::expect_equal(calls_viable, 1L) # one component
  testthat::expect_equal(calls_fit, 0L)
  testthat::expect_equal(calls_boot, 4L)   # one bootstrap per realization
})


testthat::test_that("simulate_warm component mode calls ARIMA fit helper when viable and uses returned model", {

  skip_if_not_installed("forecast")

  n <- 80L
  x <- as.numeric(stats::arima.sim(model = list(ar = 0.5), n = n))
  comps <- cbind(D1 = x)

  calls_viable <- 0L
  calls_fit <- 0L

  # Use forecast::Arima so stats::simulate dispatch works with forecast:::simulate.Arima
  fit_model <- forecast::Arima(stats::ts(x - mean(x), frequency = 1),
                               order = c(1, 0, 0),
                               include.mean = FALSE)

  testthat::local_mocked_bindings(
    .warm_arima_viable = function(x, min_n, sd_eps, min_unique) {
      calls_viable <<- calls_viable + 1L
      TRUE
    },
    .warm_fit_arima_safe = function(x, max_p, max_q, stationary, include_mean, allow_drift) {
      calls_fit <<- calls_fit + 1L
      list(model = fit_model)
    },
    .env = asNamespace("weathergenr")
  )

  out <- testthat::expect_silent(
    simulate_warm(
      components = comps,
      n = n,
      n_sim = 3L,
      seed = 42L,
      bypass_n = 25L,
      verbose = FALSE
    )
  )

  testthat::expect_equal(dim(out), c(n, 3L))
  testthat::expect_equal(calls_viable, 1L)
  testthat::expect_equal(calls_fit, 1L)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("simulate_warm constant component is carried through correctly", {

  skip_if_not_installed("forecast")

  n <- 50L
  const <- rep(10, n)
  noise <- rnorm(n)

  comps <- cbind(D1 = const, D2 = noise)

  calls_fit <- 0L

  # Again: forecast::Arima object, not stats::arima
  fit_model <- forecast::Arima(stats::ts(noise - mean(noise), frequency = 1),
                               order = c(1, 0, 0),
                               include.mean = FALSE)

  testthat::local_mocked_bindings(
    .warm_arima_viable = function(x, min_n, sd_eps, min_unique) {
      TRUE
    },
    .warm_fit_arima_safe = function(x, max_p, max_q, stationary, include_mean, allow_drift) {
      calls_fit <<- calls_fit + 1L
      list(model = fit_model)
    },
    .env = asNamespace("weathergenr")
  )

  out <- suppressWarnings(
    simulate_warm(
      components = comps,
      n = n,
      n_sim = 2L,
      seed = 7L,
      bypass_n = 25L,
      verbose = FALSE
    )
  )

  testthat::expect_equal(dim(out), c(n, 2L))
  # constant contribution should shift the mean upward
  testthat::expect_true(all(colMeans(out) > 8))
  testthat::expect_equal(calls_fit, 1L) # only D2
})

################################################################################

testthat::test_that("plot_wavelet_power and plot_wavelet_global_spectrum return plots", {
  series <- sin(seq(0, 4 * pi, length.out = 64))
  wave <- analyze_wavelet_spectrum(series, noise = "white", mode = "complete")
  signif_mask <- wave$power_signif_coi

  p <- plot_wavelet_power(
    series = series,
    period = wave$period,
    power = wave$power,
    gws = wave$gws,
    gws_signif = wave$gws_signif,
    coi = wave$coi,
    signif_mask = signif_mask,
    unit = "mm"
  )
  testthat::expect_true(inherits(p, "patchwork"))

  sim_power <- cbind(wave$gws, wave$gws * 1.1)
  p2 <- plot_wavelet_global_spectrum(
    period = wave$period,
    signif = wave$gws_signif,
    obs_power = wave$gws,
    sim_power = sim_power
  )
  testthat::expect_true(inherits(p2, "ggplot"))
})


