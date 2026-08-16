# Functions tested (relative paths):
# - R/warm_filtering.R: filter_warm_pool(), filter_warm_bounds_defaults(),
#   compute_tailmass_metrics(), compute_spectral_metrics(), relax_bounds_one_filter(),
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

testthat::test_that("compute_spectral_metrics returns cached spectra", {
  series <- sin(seq(0, 4 * pi, length.out = 32))
  sim_series_stats <- cbind(series, series * 0.9, series * 1.1)

  out <- compute_spectral_metrics(
    obs_use = series,
    sim_series_stats = sim_series_stats,
    wavelet_pars = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    ),
    cache_gws = TRUE
  )

  testthat::expect_true(all(c("period", "gws_cache", "metrics") %in% names(out)))
  testthat::expect_equal(nrow(out$gws_cache), length(out$period))
  testthat::expect_equal(ncol(out$gws_cache), ncol(sim_series_stats))
  testthat::expect_equal(length(out$metrics$spectral_cor), ncol(sim_series_stats))
  testthat::expect_true(is.matrix(out$gws_cache))
})

testthat::test_that("compute_spectral_metrics separates the display and inference signif curves", {
  # The bug: filter_warm_pool() handed the COI-masked inference curve to the
  # plot, so warm_annual_wavelet.png drew the red threshold only over the few
  # short periods the cone of influence can test and left the rest of the axis
  # bare. The display curve must be finite everywhere; the inference curve must
  # keep its NAs, because identify_significant_peaks() reads that NA as "not
  # testable" rather than "not significant".
  series <- sin(seq(0, 4 * pi, length.out = 32))
  sim_series_stats <- cbind(series, series * 0.9, series * 1.1)

  out <- compute_spectral_metrics(
    obs_use = series,
    sim_series_stats = sim_series_stats,
    wavelet_pars = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    )
  )

  testthat::expect_equal(length(out$gws_signif_display), length(out$period))
  testthat::expect_true(all(is.finite(out$gws_signif_display)))

  # Only meaningful if the two actually differ here -- otherwise the assertion
  # above would pass on a record long enough that nothing was masked.
  testthat::expect_true(anyNA(out$gws_signif))
  testthat::expect_false(identical(out$gws_signif, out$gws_signif_display))
})

testthat::test_that("plotted signif curve spans every period the spectra do", {
  # Binds the caller wiring: whatever filter_warm_pool() hands the plot must
  # render an unbroken threshold line across the full period grid.
  period <- c(2, 3, 5, 8, 12, 20)
  signif_masked <- c(0.5, 0.6, NA, NA, NA, NA)
  signif_display <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

  obs_series <- as.numeric(1:32)
  sim_series <- vapply(1:3, function(j) obs_series + j, numeric(32))
  tm <- compute_tailmass_metrics(obs_series, sim_series, 0.2, 0.8, 1e-5)

  signif_rows <- function(power_signif) {
    p <- plot_filter_diagnostics(
      obs_series = obs_series, sim_series = sim_series, pool = c(1, 2),
      rel_diff_mean = rep(0.05, 3), rel_diff_sd = rep(0.02, 3),
      tail_metrics = tm,
      power_period = period, power_obs = seq_along(period),
      power_signif = power_signif,
      wavelet_pars = list(signif_level = 0.8, noise_type = "white",
                          period_lower_limit = 2, detrend = FALSE),
      wavelet_q = c(0.5, 0.9)
    )
    # Layer 4 is the red dashed significance line.
    ggplot2::ggplot_build(p$wavelet_gws)$data[[4]]
  }

  drawn <- signif_rows(signif_display)
  testthat::expect_equal(nrow(drawn), length(period))
  testthat::expect_true(all(is.finite(drawn$y)))
  testthat::expect_equal(range(drawn$x), range(period))

  # The masked curve is what produced the truncated line; keep the contrast
  # visible so this test fails if the wiring is reverted. ggplot_build() keeps
  # the NA rows, so the truncation shows up as a missing y, not a missing row.
  masked <- suppressWarnings(signif_rows(signif_masked))
  testthat::expect_lt(max(masked$x[is.finite(masked$y)]), max(period))
})

testthat::test_that("compute_spectral_metrics keeps duplicate realizations identical", {
  series <- sin(seq(0, 6 * pi, length.out = 48))
  sim_series_stats <- cbind(series * 0.95, series * 0.95, series * 1.05)

  out <- compute_spectral_metrics(
    obs_use = series,
    sim_series_stats = sim_series_stats,
    wavelet_pars = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    ),
    cache_gws = TRUE
  )

  testthat::expect_equal(out$metrics$spectral_cor[1], out$metrics$spectral_cor[2])
  testthat::expect_equal(out$metrics$peak_match_frac[1], out$metrics$peak_match_frac[2])
  testthat::expect_equal(
    out$metrics$peak_mag_mean_abs_log_ratio[1],
    out$metrics$peak_mag_mean_abs_log_ratio[2]
  )
  testthat::expect_equal(out$gws_cache[, 1], out$gws_cache[, 2])
})

testthat::test_that("compute_spectral_metrics parallel path matches serial results", {
  set.seed(99)
  n <- 48L
  n_rlz <- 220L
  obs <- sin(seq(0, 6 * pi, length.out = n))
  sim <- replicate(n_rlz, obs + stats::rnorm(n, sd = 0.05))

  args <- list(
    obs_use = obs,
    sim_series_stats = sim,
    wavelet_pars = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    ),
    cache_gws = TRUE
  )

  serial <- do.call(compute_spectral_metrics, c(args, list(parallel = FALSE)))
  parallel_out <- do.call(compute_spectral_metrics, c(args, list(parallel = TRUE, n_cores = 2L)))

  testthat::expect_equal(serial$period, parallel_out$period)
  testthat::expect_equal(serial$gws_obs, parallel_out$gws_obs, tolerance = 1e-12)
  testthat::expect_equal(serial$gws_signif, parallel_out$gws_signif, tolerance = 1e-12)
  testthat::expect_equal(serial$gws_cache, parallel_out$gws_cache, tolerance = 1e-12)
  testthat::expect_equal(serial$metrics$spectral_cor, parallel_out$metrics$spectral_cor, tolerance = 1e-12)
  testthat::expect_equal(serial$metrics$peak_match_frac, parallel_out$metrics$peak_match_frac, tolerance = 1e-12)
  testthat::expect_equal(
    serial$metrics$peak_mag_mean_abs_log_ratio,
    parallel_out$metrics$peak_mag_mean_abs_log_ratio,
    tolerance = 1e-12
  )
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

  res_wavelet <- relax_bounds_one_filter("wavelet", b, wavelet_active_env, recompute_fn)
  testthat::expect_true(res_wavelet$changed)
})

testthat::test_that("criteria_string_compact returns filter summaries", {
  b <- filter_warm_bounds_defaults()
  tail_metrics <- list()
  spectral_diag <- list()

  out_mean <- criteria_string_compact("mean", b, tail_metrics, TRUE, spectral_diag)
  testthat::expect_true(grepl("tol", out_mean))

  out_wavelet <- criteria_string_compact("wavelet", b, tail_metrics, FALSE, spectral_diag)
  testthat::expect_identical(out_wavelet, "inactive")
})

testthat::test_that("logging helpers run without errors", {
  b <- filter_warm_bounds_defaults()
  passes <- list(mean = c(TRUE, FALSE), sd = c(TRUE, TRUE), tail_low = c(TRUE, FALSE),
                 tail_high = c(TRUE, TRUE), wavelet = c(TRUE, TRUE))
  pool <- which(Reduce("&", passes))
  tail_metrics <- compute_tailmass_metrics(c(1, 2, 3, 4), cbind(1:4, 2:5), 0.2, 0.8, 1e-5)
  spectral_diag <- list()

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
      spectral_diag = spectral_diag,
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
  testthat::expect_false("gws_cache" %in% names(out$diagnostics))
})

testthat::test_that("filter_warm_pool caches gws only when requested", {
  set.seed(2)
  obs_series <- sin(seq(0, 2 * pi, length.out = 32))
  sim_series <- sapply(1:3, function(i) obs_series + rnorm(32, sd = 0.05))

  out <- filter_warm_pool(
    obs_series = obs_series,
    sim_series = sim_series,
    n_select = 2,
    seed = 11,
    wavelet_args = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    ),
    cache_gws = TRUE,
    make_plots = FALSE,
    verbose = FALSE
  )

  testthat::expect_true("gws_cache" %in% names(out$diagnostics))
  testthat::expect_true(is.matrix(out$diagnostics$gws_cache))
  testthat::expect_equal(ncol(out$diagnostics$gws_cache), ncol(sim_series))
})

testthat::test_that("plot_filter_diagnostics returns ggplot list", {
  obs_series <- 1:32
  sim_series <- cbind(1:32, 2:33, 3:34)
  pool <- c(1, 2)

  tail_metrics <- compute_tailmass_metrics(obs_series, sim_series, 0.2, 0.8, 1e-5)
  rel_diff_mean <- c(0, 0.1, -0.1)
  rel_diff_sd <- c(0.05, -0.02, 0.01)

  power_period <- c(1, 2, 3)
  power_obs <- c(1, 2, 3)
  power_signif <- c(0.5, 0.6, 0.7)

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
    wavelet_pars = list(
      signif_level = 0.8,
      noise_type = "white",
      period_lower_limit = 2,
      detrend = FALSE
    ),
    wavelet_q = c(0.5, 0.9)
  )

  testthat::expect_true(all(c("timeseries", "stats", "wavelet_gws") %in% names(plots)))
  testthat::expect_true(inherits(plots$timeseries, "ggplot"))
  testthat::expect_true(inherits(plots$stats, "ggplot"))
  testthat::expect_true(inherits(plots$wavelet_gws, "ggplot"))

  # inherits() is close to vacuous on its own: a bare ggplot() with no data and
  # no layers passes it, and survives ggplot_build() too. Assert the structure
  # each panel is supposed to have.
  testthat::expect_equal(length(plots$timeseries$layers), 2L)
  testthat::expect_equal(length(plots$stats$layers), 2L)
  testthat::expect_equal(length(plots$wavelet_gws$layers), 4L)

  rendered <- function(p) sum(vapply(ggplot2::ggplot_build(p)$data, nrow, integer(1)))

  # The timeseries panel draws the observed record plus one line per pooled
  # realization, so its rendered row count is tied to both. This is what fails
  # if the panel stops reflecting `pool`.
  testthat::expect_equal(rendered(plots$timeseries),
                         length(obs_series) * (1L + length(pool)))

  # The stats panel draws four summary statistics for the pool and the
  # observations.
  testthat::expect_equal(rendered(plots$stats), 4L * (length(pool) + 1L))

  testthat::expect_gt(rendered(plots$wavelet_gws), 0L)
})

testthat::test_that("all three filter diagnostic panels carry the house theme", {
  # These three panels were the clearest case of the drift: two used
  # theme_light() and the stats panel used theme_bw(), inside one function
  # returning them as a set. panel.border$colour discriminates all three themes
  # the package used to mix; panel.background$fill would not, since theme_light
  # and theme_bw both give white.
  obs_series <- as.numeric(1:32)
  sim_series <- vapply(1:3, function(j) obs_series + j, numeric(32))
  tm <- compute_tailmass_metrics(obs_series, sim_series, 0.2, 0.8, 1e-5)

  plots <- plot_filter_diagnostics(
    obs_series = obs_series, sim_series = sim_series, pool = c(1, 2),
    rel_diff_mean = rep(0.05, 3), rel_diff_sd = rep(0.02, 3),
    tail_metrics = tm,
    power_period = c(1, 2, 3), power_obs = c(1, 2, 3),
    power_signif = c(0.5, 0.6, 0.7),
    wavelet_pars = list(signif_level = 0.8, noise_type = "white",
                        period_lower_limit = 2, detrend = FALSE),
    wavelet_q = c(0.5, 0.9)
  )

  house <- theme_weathergenr()
  for (nm in c("timeseries", "stats", "wavelet_gws")) {
    testthat::expect_identical(plots[[nm]]$theme$panel.border$colour,
                               house$panel.border$colour, info = nm)
    testthat::expect_identical(plots[[nm]]$theme$text$size,
                               house$text$size, info = nm)
  }
})

testthat::test_that("plot_filter_diagnostics tracks the pool it is given", {
  # The row-count identities above only bind if they move with their inputs;
  # a panel that ignored `pool` would satisfy a single fixed number.
  obs_series <- as.numeric(1:32)
  sim_series <- vapply(1:5, function(j) obs_series + j, numeric(32))

  build <- function(pool) {
    tm <- compute_tailmass_metrics(obs_series, sim_series, 0.2, 0.8, 1e-5)
    p <- plot_filter_diagnostics(
      obs_series = obs_series, sim_series = sim_series, pool = pool,
      rel_diff_mean = rep(0.05, 5), rel_diff_sd = rep(0.02, 5),
      tail_metrics = tm,
      power_period = c(1, 2, 3), power_obs = c(1, 2, 3),
      power_signif = c(0.5, 0.6, 0.7),
      wavelet_pars = list(signif_level = 0.8, noise_type = "white",
                          period_lower_limit = 2, detrend = FALSE),
      wavelet_q = c(0.5, 0.9)
    )
    sum(vapply(ggplot2::ggplot_build(p$timeseries)$data, nrow, integer(1)))
  }

  testthat::expect_equal(build(c(1, 2)), 32L * 3L)
  testthat::expect_equal(build(c(1, 2, 3)), 32L * 4L)
  testthat::expect_equal(build(1:5), 32L * 6L)
})


################################################################################
# Fork equivalence: .compute_gws_batch() vs analyze_wavelet_spectrum()
################################################################################

testthat::test_that(".compute_gws_batch reproduces analyze_wavelet_spectrum's GWS", {
  # .compute_gws_batch() is a batched reimplementation of the CWT in
  # analyze_wavelet_spectrum(), kept separate for speed over large ensembles.
  # Nothing but this test holds the two together: a change to dj, s0, k0, the
  # zero-padding, the detrend or the variance normalisation in one must be
  # mirrored in the other.
  set.seed(3)
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), 40)) * 20 + 800

  for (dt in c(FALSE, TRUE)) {
    ref <- analyze_wavelet_spectrum(x, signif = 0.90, noise = "red",
                                    min_period = 2, detrend = dt, mode = "complete")
    fork <- weathergenr:::.compute_gws_batch(
      sim_matrix   = matrix(x, ncol = 1),
      wavelet_pars = list(detrend = dt),
      period       = NULL
    )

    testthat::expect_equal(nrow(fork), length(ref$gws_unmasked))
    testthat::expect_equal(as.numeric(fork[, 1]), as.numeric(ref$gws_unmasked),
                           tolerance = 1e-10)
  }
})

testthat::test_that(".compute_gws_batch is column-wise independent", {
  # The batching must not let one trace influence another: column j of the
  # result has to equal the single-series result for that trace.
  set.seed(4)
  m <- cbind(
    as.numeric(stats::arima.sim(list(ar = 0.3), 40)) * 10 + 500,
    as.numeric(stats::arima.sim(list(ar = 0.7), 40)) * 30 + 900
  )

  batch <- weathergenr:::.compute_gws_batch(m, list(detrend = TRUE), period = NULL)

  for (j in 1:2) {
    one <- weathergenr:::.compute_gws_batch(m[, j, drop = FALSE],
                                           list(detrend = TRUE), period = NULL)
    testthat::expect_equal(as.numeric(batch[, j]), as.numeric(one[, 1]),
                           tolerance = 1e-12)
  }

  # Chunking must not change the answer either.
  chunked <- weathergenr:::.compute_gws_batch(m, list(detrend = TRUE),
                                              period = NULL, chunk_size = 1L)
  testthat::expect_equal(chunked, batch, tolerance = 1e-12)
})


################################################################################
# Deterministic length harmonisation
################################################################################

testthat::test_that("the observed benchmark window does not depend on the seed", {
  set.seed(1)
  obs <- as.numeric(stats::arima.sim(list(ar = 0.5), 40)) * 30 + 900
  sim <- matrix(stats::rnorm(25 * 60, mean = 900, sd = 30), nrow = 25)

  runs <- lapply(c(1L, 2L, 999L), function(s) {
    filter_warm_pool(obs_series = obs, sim_series = sim, n_select = 3L,
                     seed = s, make_plots = FALSE, verbose = FALSE)
  })

  # The window is the tail of the observed record, identically for every seed.
  for (r in runs) {
    testthat::expect_equal(unname(r$diagnostics$obs_window),
                           c(length(obs) - 25L + 1L, length(obs)))
  }

  # The acceptance statistics derive from that window, so they must agree too.
  # (Which traces are drawn from the pool still depends on the seed; that is
  # ensemble sampling, not a moving criterion.)
  pools <- vapply(runs, function(r) length(r$diagnostics$pool_idx), integer(1))
  testthat::expect_equal(length(unique(pools)), 1L)
})

testthat::test_that("every candidate is scored on the same simulated window", {
  # n_sim > n_obs is the direction the packaged baseline never exercises: each
  # candidate used to get its own random window, so candidates were not
  # compared on the same footing.
  set.seed(2)
  obs <- as.numeric(stats::arima.sim(list(ar = 0.5), 20)) * 30 + 900
  sim <- matrix(stats::rnorm(60 * 40, mean = 900, sd = 30), nrow = 60)

  res <- filter_warm_pool(obs_series = obs, sim_series = sim, n_select = 3L,
                          seed = 7L, make_plots = FALSE, verbose = FALSE)

  wi <- res$diagnostics$window_index
  testthat::expect_equal(length(wi), 40L)
  testthat::expect_equal(length(unique(lapply(wi, unname))), 1L)
  testthat::expect_equal(unname(wi[[1]]), c(60L - 20L + 1L, 60L))

  # Selected traces are returned at full length, not windowed.
  testthat::expect_equal(nrow(res$selected), 60L)
})

testthat::test_that("compute_tailmass_metrics normalises each side by its own length", {
  set.seed(3)
  obs <- stats::rnorm(20, mean = 900, sd = 30)

  # The same trace repeated end-to-end has the same per-observation tail mass,
  # so doubling its length must not double the metric.
  one <- matrix(stats::rnorm(20, mean = 900, sd = 30), ncol = 1)
  two <- matrix(rep(one[, 1], 2), ncol = 1)

  m1 <- compute_tailmass_metrics(obs, one, tail_low_p = 0.2,
                                 tail_high_p = 0.8, tail_eps = 1e-5)
  m2 <- compute_tailmass_metrics(obs, two, tail_low_p = 0.2,
                                 tail_high_p = 0.8, tail_eps = 1e-5)

  testthat::expect_equal(m2$M_sim_low, m1$M_sim_low)
  testthat::expect_equal(m2$M_sim_high, m1$M_sim_high)

  # The observed side is unaffected by the simulated length.
  testthat::expect_equal(m2$M_obs_low, m1$M_obs_low)
})
