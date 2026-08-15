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

testthat::test_that("internal ARMA matrix simulation matches full-matrix reference", {
  ref_simulate_arma_matrix <- function(ar, ma, sd, n, n_sim, burnin = 100L) {
    ar <- as.numeric(ar)
    ma <- as.numeric(ma)
    sd <- as.numeric(sd)
    n <- as.integer(n)
    n_sim <- as.integer(n_sim)
    burnin <- as.integer(burnin)

    p <- length(ar)
    q <- length(ma)
    n_tot <- n + burnin
    max_lag <- max(p, q, 1L)

    eps <- matrix(
      stats::rnorm((n_tot + max_lag) * n_sim, mean = 0, sd = sd),
      nrow = n_tot + max_lag,
      ncol = n_sim
    )
    x <- matrix(0, nrow = n_tot + max_lag, ncol = n_sim)

    for (t in (max_lag + 1L):(n_tot + max_lag)) {
      xt <- eps[t, ]
      if (q >= 1L) xt <- xt + ma[1L] * eps[t - 1L, ]
      if (q >= 2L) xt <- xt + ma[2L] * eps[t - 2L, ]
      if (p >= 1L) xt <- xt + ar[1L] * x[t - 1L, ]
      if (p >= 2L) xt <- xt + ar[2L] * x[t - 2L, ]
      if (p > 2L) {
        for (i in 3L:p) xt <- xt + ar[i] * x[t - i, ]
      }
      if (q > 2L) {
        for (i in 3L:q) xt <- xt + ma[i] * eps[t - i, ]
      }
      x[t, ] <- xt
    }

    x[(max_lag + burnin + 1L):(max_lag + burnin + n), , drop = FALSE]
  }

  pars <- list(ar = c(0.6, -0.2), ma = c(0.3), sd = 1.1, n = 24L, n_sim = 1005L, burnin = 20L)

  set.seed(123)
  out_new <- do.call(weathergenr:::.simulate_arma_matrix, pars)
  set.seed(123)
  out_ref <- do.call(ref_simulate_arma_matrix, pars)

  testthat::expect_equal(out_new, out_ref, tolerance = 0)
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
    .default_block_len = function(n) 5L
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

test_that("simulate_warm bypass mode caps ARMA order for very short fits", {
  recorded_max_p <- integer(0)
  recorded_max_q <- integer(0)

  testthat::local_mocked_bindings(
    .fit_warm_arima_forecast = function(x, max_p, max_q) {
      recorded_max_p <<- c(recorded_max_p, as.integer(max_p))
      recorded_max_q <<- c(recorded_max_q, as.integer(max_q))
      NULL
    },
    .block_bootstrap = function(x, n, block_len) rep(mean(x), n),
    .default_block_len = function(n) 5L
  )

  simulate_warm(
    n = 14L,
    n_sim = 2L,
    series_obs = rnorm(14),
    bypass_n = 15L,
    verbose = FALSE
  )
  simulate_warm(
    n = 20L,
    n_sim = 2L,
    series_obs = rnorm(20),
    bypass_n = 25L,
    verbose = FALSE
  )

  expect_equal(recorded_max_p, c(1L, 2L))
  expect_equal(recorded_max_q, c(1L, 2L))
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
    }
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
    }
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
    }
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

test_that("simulate_warm warns on large pre-correction variance mismatch with diagnostics", {
  n <- 30L
  comp1 <- seq_len(n)
  comp2 <- 0.5 * seq_len(n)
  components <- cbind(D1 = comp1, D2 = comp2)

  testthat::local_mocked_bindings(
    .warm_arima_viable = function(x, min_n, sd_eps, min_unique) TRUE,
    .warm_fit_arima_safe = function(x, max_p, max_q, stationary, include_mean, allow_drift) {
      list(model = "mock_fit", residuals = NULL)
    },
    .warm_simulate_from_fit = function(fit, n, n_sim) matrix(0, nrow = n, ncol = n_sim)
  )

  expect_warning(
    simulate_warm(
      components = components,
      n = n,
      n_sim = 2L,
      seed = 1L,
      check_diagnostics = TRUE,
      verbose = FALSE
    ),
    "pre-correction variance mismatch is high"
  )
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



# ==============================================================================
# simulate_warm(): argument validation
#
# Every check is an early exit, so these are effectively free. They run in
# declaration order, which is why each case below supplies valid values for the
# arguments validated ahead of the one under test.
# ==============================================================================

testthat::test_that("simulate_warm validates its arguments", {
  ok <- as.numeric(1:40)

  testthat::expect_error(
    simulate_warm(components = NULL, n = NULL, series_obs = ok),
    "Input 'n' must be specified"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 0L, series_obs = ok),
    "'n' must be a positive integer"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 0L, series_obs = ok),
    "'n_sim' must be a positive integer"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, bypass_n = 4L, series_obs = ok),
    "'bypass_n' must be an integer >= 5"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, verbose = c(TRUE, FALSE), series_obs = ok),
    "'verbose' must be TRUE or FALSE"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, match_variance = "yes", series_obs = ok),
    "'match_variance' must be TRUE or FALSE"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, var_tol = 1.5, series_obs = ok),
    "'var_tol' must be between 0 and 1"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, check_diagnostics = 1, series_obs = ok),
    "'check_diagnostics' must be TRUE or FALSE"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, seed = "abc", series_obs = ok),
    "'seed' must be NULL or a single finite number"
  )
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, series_obs = c(1, NA, 3)),
    "'series_obs' must be numeric with no missing values"
  )

  # neither series_obs nor components: the fit length cannot be established
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L),
    "Cannot determine observed series length"
  )

  # component mode reached with NULL components is rejected separately
  testthat::expect_error(
    simulate_warm(components = NULL, n = 40L, n_sim = 2L, series_obs = ok, bypass_n = 5L),
    "Input 'components' must not be NULL for component mode"
  )
})

testthat::test_that("simulate_warm rejects malformed components in component mode", {
  n <- 40L
  good <- as.numeric(1:n)

  testthat::expect_error(
    simulate_warm(components = matrix("a", nrow = n, ncol = 2), n = n, n_sim = 2L,
                  series_obs = good, verbose = FALSE),
    "Matrix 'components' must be numeric"
  )
  testthat::expect_error(
    simulate_warm(components = matrix(1, nrow = n - 1L, ncol = 2), n = n, n_sim = 2L,
                  series_obs = good, verbose = FALSE),
    "Matrix 'components' must have n rows"
  )
  testthat::expect_error(
    simulate_warm(components = data.frame(a = as.character(seq_len(n)), b = 1),
                  n = n, n_sim = 2L, series_obs = good, verbose = FALSE),
    "All columns of 'components' must be numeric"
  )
  testthat::expect_error(
    simulate_warm(components = list(a = as.character(seq_len(n))),
                  n = n, n_sim = 2L, series_obs = good, verbose = FALSE),
    "All elements of 'components' must be numeric vectors"
  )
  testthat::expect_error(
    simulate_warm(components = list(a = good, b = good[-1]),
                  n = n, n_sim = 2L, series_obs = good, verbose = FALSE),
    "All component vectors must have length 'n'"
  )
  testthat::expect_error(
    simulate_warm(components = "not-components", n = n, n_sim = 2L,
                  series_obs = good, verbose = FALSE),
    "'components' must be a matrix, data.frame, or list"
  )

  na_comp <- cbind(good, c(NA_real_, good[-1]))
  testthat::expect_error(
    simulate_warm(components = na_comp, n = n, n_sim = 2L, series_obs = good, verbose = FALSE),
    "Missing values detected in component"
  )
})

# ==============================================================================
# simulate_warm(): bypass mode, ARMA-success path
#
# The existing bypass tests above drive the block-bootstrap fallback (fit
# fails). These cover the other branch: a short series with enough structure
# for .fit_warm_arima_forecast() to return a model, which is the only route
# through .warm_simulate_from_fit() and the variance-matching step.
# ==============================================================================

testthat::test_that("simulate_warm bypass mode simulates from a successful ARMA fit", {
  set.seed(11)
  # 12 points (< default bypass_n = 15) with clear AR(1) structure
  series <- as.numeric(stats::filter(rnorm(12, sd = 1), 0.6, method = "recursive")) + 20

  sim <- simulate_warm(
    components = NULL,
    n = 12L,
    n_sim = 25L,
    series_obs = series,
    seed = 42,
    verbose = FALSE
  )

  testthat::expect_true(is.matrix(sim))
  testthat::expect_equal(dim(sim), c(12L, 25L))
  testthat::expect_true(all(is.finite(sim)))

  # the simulation is centred on the observed mean, not on zero
  testthat::expect_lt(abs(mean(sim) - mean(series)), 3 * stats::sd(series))

  # and it is reproducible for a fixed seed
  sim2 <- simulate_warm(
    components = NULL, n = 12L, n_sim = 25L,
    series_obs = series, seed = 42, verbose = FALSE
  )
  testthat::expect_identical(sim, sim2)
})

testthat::test_that("simulate_warm bypass mode reconstructs the series from components", {
  set.seed(12)
  n <- 12L
  comp1 <- as.numeric(stats::filter(rnorm(n, sd = 1), 0.5, method = "recursive"))
  comp2 <- seq(0, 2, length.out = n)
  comps <- cbind(comp1, comp2)

  # With no series_obs, bypass mode must rebuild it as the component rowSums.
  from_comps <- simulate_warm(
    components = comps, n = n, n_sim = 10L, seed = 5, verbose = FALSE
  )
  from_series <- simulate_warm(
    components = NULL, n = n, n_sim = 10L,
    series_obs = rowSums(comps), seed = 5, verbose = FALSE
  )

  testthat::expect_equal(dim(from_comps), c(n, 10L))
  testthat::expect_identical(from_comps, from_series)
})

testthat::test_that("simulate_warm bypass mode reports itself when verbose", {
  set.seed(13)
  series <- as.numeric(stats::filter(rnorm(12, sd = 1), 0.6, method = "recursive"))

  testthat::expect_message(
    simulate_warm(components = NULL, n = 12L, n_sim = 5L,
                  series_obs = series, seed = 1, verbose = TRUE),
    "Bypass mode"
  )
})

# ==============================================================================
# Internal helpers reached only from simulate_warm()
# ==============================================================================

testthat::test_that(".reconstruct_series_from_components sums across every supported shape", {
  n <- 6L
  a <- as.numeric(1:n)
  b <- as.numeric((1:n) * 10)
  expected <- a + b

  testthat::expect_identical(
    weathergenr:::.reconstruct_series_from_components(cbind(a, b), n = n),
    expected
  )
  testthat::expect_identical(
    weathergenr:::.reconstruct_series_from_components(data.frame(a = a, b = b), n = n),
    expected
  )
  testthat::expect_identical(
    weathergenr:::.reconstruct_series_from_components(list(a = a, b = b), n = n),
    expected
  )

  testthat::expect_error(
    weathergenr:::.reconstruct_series_from_components(cbind(a, b), n = n + 1L),
    "'components' must have n rows"
  )
  testthat::expect_error(
    weathergenr:::.reconstruct_series_from_components(list(a = a, b = b[-1]), n = n),
    "all component vectors must have length n"
  )
  testthat::expect_error(
    weathergenr:::.reconstruct_series_from_components("nope", n = n),
    "unsupported 'components' type"
  )
})

testthat::test_that(".warm_arima_viable screens series that cannot support an ARMA fit", {
  viable <- as.numeric(stats::filter(rnorm(40), 0.5, method = "recursive"))

  testthat::expect_true(weathergenr:::.warm_arima_viable(viable))

  testthat::expect_false(weathergenr:::.warm_arima_viable("a"))
  testthat::expect_false(weathergenr:::.warm_arima_viable(c(1, NA, 3, 4, 5, 6, 7, 8, 9)))
  testthat::expect_false(weathergenr:::.warm_arima_viable(as.numeric(1:5)))   # shorter than min_n
  testthat::expect_false(weathergenr:::.warm_arima_viable(rep(1, 20)))        # zero variance
  testthat::expect_false(
    weathergenr:::.warm_arima_viable(rep(c(1, 2), each = 10))                 # too few unique
  )
})

testthat::test_that(".fit_warm_arima_forecast delegates to the safe ARMA fitter", {
  set.seed(14)
  x <- as.numeric(stats::filter(rnorm(60), 0.6, method = "recursive"))
  x <- x - mean(x)

  fit <- weathergenr:::.fit_warm_arima_forecast(x, max_p = 2L, max_q = 2L)

  testthat::expect_type(fit, "list")
  testthat::expect_true(all(c("ar", "ma", "sigma2", "residuals", "order") %in% names(fit)))
  testthat::expect_true(is.finite(fit$sigma2))

  # the wrapper adds nothing of its own beyond forwarding its arguments
  testthat::expect_identical(
    fit,
    weathergenr:::.warm_fit_arima_safe(x, max_p = 2L, max_q = 2L)
  )

  # an unfittable series propagates the NULL rather than erroring
  testthat::expect_null(weathergenr:::.fit_warm_arima_forecast(c(1, 2), max_p = 1L, max_q = 1L))
})

testthat::test_that("simulate_warm accepts data.frame and list components in component mode", {
  set.seed(15)
  n <- 40L
  comp1 <- as.numeric(stats::filter(rnorm(n), 0.5, method = "recursive"))
  comp2 <- as.numeric(stats::filter(rnorm(n), 0.2, method = "recursive"))

  as_matrix <- simulate_warm(
    components = cbind(comp1, comp2), n = n, n_sim = 5L, seed = 3, verbose = FALSE
  )
  as_df <- simulate_warm(
    components = data.frame(comp1 = comp1, comp2 = comp2),
    n = n, n_sim = 5L, seed = 3, verbose = FALSE
  )
  as_list <- simulate_warm(
    components = list(comp1 = comp1, comp2 = comp2),
    n = n, n_sim = 5L, seed = 3, verbose = FALSE
  )

  testthat::expect_equal(dim(as_matrix), c(n, 5L))
  testthat::expect_equal(as_df, as_matrix)
  testthat::expect_equal(as_list, as_matrix)
})


################################################################################
# .fast_col_sd / .variance_match_matrix
################################################################################

testthat::test_that(".fast_col_sd returns the sample standard deviation", {
  set.seed(1)
  m <- matrix(stats::rnorm(40 * 5), nrow = 40)

  # Must match stats::sd, because the targets it is compared against and
  # rescaled to are produced by stats::sd. A population sd here left every
  # corrected column at target_sd * sqrt(n / (n - 1)).
  expect_equal(weathergenr:::.fast_col_sd(m), apply(m, 2, stats::sd))

  # Fewer than two rows cannot support an estimate.
  expect_true(all(is.na(weathergenr:::.fast_col_sd(matrix(1, nrow = 1, ncol = 3)))))
})

testthat::test_that(".variance_match_matrix rescales onto the target without bias", {
  set.seed(2)
  n <- 40
  s <- matrix(stats::rnorm(n * 200, mean = 900, sd = 40), nrow = n)
  target <- 69.371

  out <- weathergenr:::.variance_match_matrix(s, target_sd = target, tol = 0.1)
  sds <- apply(out, 2, stats::sd)

  # Every column here is well outside tol, so all are corrected and must land
  # exactly on the target -- not on target * sqrt(n / (n - 1)).
  expect_equal(unname(sds), rep(target, ncol(s)))
  expect_false(any(abs(sds - target * sqrt(n / (n - 1))) < 1e-8))
})

testthat::test_that(".variance_match_matrix preserves each column's own mean", {
  set.seed(3)
  n <- 40
  s <- matrix(stats::rnorm(n * 200, mean = 900, sd = 40), nrow = n)

  out <- weathergenr:::.variance_match_matrix(s, target_sd = 69.371, tol = 0.1)

  # Rescaling is about the column's own mean. Re-centring on a scalar target
  # put a point mass on it and made filter_warm_pool()'s mean criterion unable
  # to reject any corrected realization.
  expect_equal(colMeans(out), colMeans(s))
  expect_gt(length(unique(round(colMeans(out), 9))), 190)
})

testthat::test_that("simulate_warm keeps spread in the trace mean in component mode", {
  # A single variable component keeps nvar < 2, so Cholesky decorrelation is
  # skipped and the per-component variance-matching branch runs -- the branch
  # the packaged fixture never reaches.
  set.seed(11)
  comps <- cbind(
    D1 = as.numeric(stats::arima.sim(list(ar = 0.5), 40)) * 30 + 400,
    S1 = rep(500, 40)
  )
  obs_total <- rowSums(comps)

  sim <- simulate_warm(components = comps, n = 40, n_sim = 300, seed = 5,
                       match_variance = TRUE, verbose = FALSE)

  expect_true(all(is.finite(sim)))

  # The per-component means must be added back exactly once.
  expect_equal(mean(sim), mean(obs_total), tolerance = 0.01)

  # No realization's mean is pinned to a shared value.
  expect_gt(length(unique(round(colMeans(sim), 9))), 290)
  expect_gt(stats::sd(colMeans(sim)), 0)
})
