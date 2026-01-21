# Functions tested (relative paths):
# - R/warm_filtering.R: filter_warm_pool(), filter_warm_bounds_defaults(),
#   compute_tailmass_metrics(), compute_wavelet_metrics(), relax_bounds_one_filter(),
#   criteria_string_compact(), log_filtering_start(), log_filter_iteration(),
#   log_final_summary()
# - R/warm_filtering_plots.R: plot_filter_diagnostics()

testthat::test_that("filter_warm_bounds_defaults returns expected keys", {
  b <- filter_warm_bounds_defaults()
  testthat::expect_true(all(c("mean", "sd", "tail_low_p", "tail_high_p",
                              "tail_tol_log", "tail_eps", "relax_max_iter") %in% names(b)))
  testthat::expect_true(is.numeric(b$mean))
  testthat::expect_true(is.numeric(b$sd))
  testthat::expect_true(is.numeric(b$tail_low_p))
  testthat::expect_true(is.numeric(b$tail_high_p))
})

testthat::test_that("compute_tailmass_metrics returns finite vectors", {
  obs_use <- c(1, 2, 3, 4)
  sim_series_stats <- cbind(c(1, 2, 2, 3), c(2, 3, 4, 5))

  out <- compute_tailmass_metrics(
    obs_use = obs_use,
    sim_series_stats = sim_series_stats,
    tail_low_p = 0.2,
    tail_high_p = 0.8,
    tail_eps = 1e-5
  )

  testthat::expect_true(is.finite(out$thr_low))
  testthat::expect_true(is.finite(out$thr_high))
  testthat::expect_equal(length(out$M_sim_low), ncol(sim_series_stats))
  testthat::expect_equal(length(out$logdiff_high), ncol(sim_series_stats))
})

testthat::test_that("compute_wavelet_metrics returns cached spectra", {
  series <- sin(seq(0, 4 * pi, length.out = 32))
  sim_series_stats <- cbind(series, series * 0.9, series * 1.1)

  out <- compute_wavelet_metrics(
    obs_use = series,
    sim_series_stats = sim_series_stats,
    wavelet_pars = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    ),
    padding = TRUE,
    min_bg = 1e-12
  )

  testthat::expect_true(all(c("power_period", "gws_cache", "gws_cache_unmasked") %in% names(out)))
  testthat::expect_equal(nrow(out$gws_cache), length(out$power_period))
  testthat::expect_equal(ncol(out$gws_cache), ncol(sim_series_stats))
})

testthat::test_that("relax_bounds_one_filter updates bounds as expected", {
  b_list <- filter_warm_bounds_defaults()
  b <- list2env(b_list, parent = environment())

  wavelet_active_env <- new.env(parent = emptyenv())
  assign("wavelet_active", TRUE, envir = wavelet_active_env)

  recompute_called <- FALSE
  recompute_fn <- function() recompute_called <<- TRUE

  res_mean <- relax_bounds_one_filter("mean", b, wavelet_active_env, recompute_fn)
  testthat::expect_true(res_mean$changed)

  b$tail_tol_log <- b$relax_tail_tol_log_max
  b$tail_low_p <- min(b$tail_low_p, b$relax_tail_p_low_max - 0.05)
  res_tail <- relax_bounds_one_filter("tail_low", b, wavelet_active_env, recompute_fn)
  testthat::expect_true(res_tail$changed)
  testthat::expect_true(recompute_called)

  b$sig_frac <- b$relax_wavelet_sig_frac_min + 0.10
  res_wavelet <- relax_bounds_one_filter("wavelet", b, wavelet_active_env, recompute_fn)
  testthat::expect_true(res_wavelet$changed)
})

testthat::test_that("criteria_string_compact returns filter summaries", {
  b <- filter_warm_bounds_defaults()
  tail_metrics <- list()
  wavelet_pars <- list()

  out_mean <- criteria_string_compact("mean", b, tail_metrics, TRUE, wavelet_pars)
  testthat::expect_true(grepl("tol", out_mean))

  out_wavelet <- criteria_string_compact("wavelet", b, tail_metrics, FALSE, wavelet_pars)
  testthat::expect_identical(out_wavelet, "inactive")
})

testthat::test_that("logging helpers run without errors", {
  b <- filter_warm_bounds_defaults()
  passes <- list(mean = c(TRUE, FALSE), sd = c(TRUE, TRUE), tail_low = c(TRUE, FALSE),
                 tail_high = c(TRUE, TRUE), wavelet = c(TRUE, TRUE))
  pool <- which(Reduce("&", passes))
  tail_metrics <- compute_tailmass_metrics(c(1, 2, 3, 4), cbind(1:4, 2:5), 0.2, 0.8, 1e-5)
  wavelet_pars <- list(signif_level = 0.8, noise_type = "white", period_lower_limit = 2, detrend = FALSE)

  testthat::expect_silent(suppressMessages(
    log_filtering_start(4, 4, 2, 1, c("mean", "sd", "tail_low", "tail_high", "wavelet"))
  ))

  testthat::expect_silent(suppressMessages(
    log_filter_iteration(
      iter = 0L,
      passes = passes,
      pool = pool,
      n_total = 2,
      target = 1,
      bounds = b,
      tail_metrics = tail_metrics,
      wavelet_active = TRUE,
      wavelet_pars = wavelet_pars,
      note = "test"
    )
  ))

  testthat::expect_silent(suppressMessages(
    log_final_summary(pool_size = 1, n_total = 2, n_sampled = 1, relaxation_level = "test")
  ))
})

testthat::test_that("filter_warm_pool returns expected structure", {
  set.seed(1)
  obs_series <- sin(seq(0, 2 * pi, length.out = 32))
  sim_series <- sapply(1:4, function(i) obs_series + rnorm(32, sd = 0.1))

  out <- filter_warm_pool(
    obs_series = obs_series,
    sim_series = sim_series,
    n_select = 2,
    seed = 10,
    wavelet_args = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    ),
    make_plots = FALSE,
    verbose = FALSE
  )

  testthat::expect_true(all(c("pool", "selected", "summary", "diagnostics", "plots") %in% names(out)))
  testthat::expect_equal(ncol(out$selected), 2)
  testthat::expect_equal(nrow(out$summary), 5)
  testthat::expect_equal(length(out$diagnostics$selected_idx), 2)
})

testthat::test_that("plot_filter_diagnostics returns ggplot list", {
  obs_series <- 1:10
  sim_series <- cbind(1:10, 2:11, 3:12)
  pool <- c(1, 2)

  tail_metrics <- compute_tailmass_metrics(obs_series, sim_series, 0.2, 0.8, 1e-5)
  rel_diff_mean <- c(0, 0.1, -0.1)
  rel_diff_sd <- c(0.05, -0.02, 0.01)

  power_period <- c(1, 2, 3)
  power_obs <- c(1, 2, 3)
  power_signif <- c(0.5, 0.6, 0.7)
  gws_cache <- matrix(c(1, 2, 3, 1.1, 2.1, 3.1, 0.9, 1.9, 2.9), nrow = 3)

  plots <- plot_filter_diagnostics(
    obs_series = obs_series,
    sim_series = sim_series,
    pool = pool,
    rel_diff_mean = rel_diff_mean,
    rel_diff_sd = rel_diff_sd,
    tail_metrics = tail_metrics,
    power_period = power_period,
    power_obs = power_obs,
    power_signif = power_signif,
    gws_cache = gws_cache,
    wavelet_q = c(0.5, 0.9)
  )

  testthat::expect_true(all(c("timeseries", "stats", "wavelet_gws") %in% names(plots)))
  testthat::expect_true(inherits(plots$timeseries, "ggplot"))
  testthat::expect_true(inherits(plots$stats, "ggplot"))
  testthat::expect_true(inherits(plots$wavelet_gws, "ggplot"))
})
