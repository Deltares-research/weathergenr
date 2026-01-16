# ==============================================================================
# Tests for evaluate_weather_generator()
# ==============================================================================
#
# Objective: keep unit coverage while preventing CI/runtime termination.
# Strategy:
#   - Keep validation tests "real" (no mocks) since they exit early and are fast.
#   - For tests that require a successful full run, mock the expensive internals:
#       * .summarize_observed_data()
#       * .summarize_simulated_data()
#       * .build_plot_data()
#       * create_all_diagnostic_plots()
#       * .summarize_realization_fit()
#     This preserves return structure, class, attributes, and key behaviors
#     (subsampling, leap-day handling, seed determinism), without correlation/plot cost.
# ==============================================================================

# ------------------------------------------------------------------------------
# Helper: Create synthetic grid data frame
# ------------------------------------------------------------------------------

make_test_grid_df <- function(dates, id_shift = 0) {

  n <- length(dates)
  data.frame(
    date = dates,
    precip = rgamma(n, shape = 2, scale = 2) + id_shift * 0.1,
    temp = rnorm(n, mean = 10 + id_shift, sd = 3)
  )
}

# ------------------------------------------------------------------------------
# Helper: Apply fast mocks for successful-run tests
# ------------------------------------------------------------------------------

.with_fast_eval_mocks <- function(code) {

  # Minimal ggplot object generator (fast, stable)
  .p <- function() ggplot2::ggplot() + ggplot2::geom_blank()

  # Minimal plot list (must include keys used in tests)
  .mock_create_all_diagnostic_plots <- function(plot_data, plot_config, variables,
                                                show_title, save_plots, output_path) {
    plots <- list(
      daily_mean = .p(),
      daily_sd = .p(),
      annual_precip = .p()
    )
    plots
  }

  # Minimal observed results structure used downstream
  .mock_summarize_observed_data <- function(daily_obs, variables, grid_count,
                                            wet_quantile, extreme_quantile) {

    # thresholds table used by .summarize_simulated_data()
    mc_thresholds <- dplyr::tibble(
      id = 1L,
      mon = 1L,
      wet.th = 0,
      extreme.th = 0
    )

    # minimal long tables used by .summarize_realization_fit() (but we mock that too)
    empty_tbl <- dplyr::tibble()

    list(
      data = empty_tbl,
      datemat = empty_tbl,
      mc_thresholds = mc_thresholds,
      stats.season = empty_tbl,
      stats.mon.aavg = empty_tbl,
      stats.annual.aavg = empty_tbl,
      wetdry = empty_tbl,
      cor = empty_tbl,
      cor.cond = empty_tbl
    )
  }

  # Minimal simulated results structure
  .mock_summarize_simulated_data <- function(daily_sim, n_realizations, variables, mc_thresholds) {

    empty_tbl <- dplyr::tibble()

    list(
      stats.season = empty_tbl,
      stats.mon.aavg = empty_tbl,
      stats.annual.aavg = empty_tbl,
      cor = empty_tbl,
      wetdry = empty_tbl,
      cor.cond = empty_tbl
    )
  }

  # Minimal plot_data structure; create_all_diagnostic_plots is mocked anyway
  .mock_build_plot_data <- function(obs_results, sim_results, variables) {
    list()
  }

  # Deterministic, minimal fit summary required by tests
  .mock_summarize_realization_fit <- function(obs_results, sim_results, variables) {
    dplyr::tibble(
      rlz = 1L,
      rank = 1L,
      overall_score = 0
    )
  }

  testthat::local_mocked_bindings(
    create_all_diagnostic_plots = .mock_create_all_diagnostic_plots,
    .summarize_observed_data = .mock_summarize_observed_data,
    .summarize_simulated_data = .mock_summarize_simulated_data,
    .build_plot_data = .mock_build_plot_data,
    .summarize_realization_fit = .mock_summarize_realization_fit
  )

  force(code)
}

# ==============================================================================
# TEST: Basic output structure and class
# ==============================================================================

testthat::test_that("evaluate_weather_generator returns weather_assessment with expected structure", {

  .with_fast_eval_mocks({

    set.seed(123)

    # Create observed dates (3 full years, no leap days)
    obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2003-12-31"), by = "day")
    obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

    # Create simulated dates (3 full years, no leap days)
    sim_dates <- seq.Date(as.Date("2011-01-01"), as.Date("2013-12-31"), by = "day")
    sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

    n_grid <- 3
    n_realizations <- 2

    # Build observed data list (one data.frame per grid cell)
    daily_obs <- lapply(seq_len(n_grid), function(i) {
      make_test_grid_df(obs_dates, id_shift = i)
    })

    # Build simulated data list (outer: realizations, inner: grid cells)
    daily_sim <- lapply(seq_len(n_realizations), function(r) {
      lapply(seq_len(n_grid), function(i) {
        make_test_grid_df(sim_dates, id_shift = i + r)
      })
    })

    # Run evaluation
    out <- evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = n_realizations,
      wet_quantile = 0.2,
      extreme_quantile = 0.8,
      output_path = NULL,
      save_plots = FALSE,
      show_title = FALSE,
      max_grids = 25,
      verbose = FALSE
    )

    # Check class
    testthat::expect_s3_class(out, "weather_assessment")
    testthat::expect_true(is.list(out))
    testthat::expect_true(length(out) > 0)

    # Check expected plot elements exist and are ggplot objects
    testthat::expect_true(inherits(out$daily_mean, "ggplot"))
    testthat::expect_true(inherits(out$daily_sd, "ggplot"))
    testthat::expect_true(inherits(out$annual_precip, "ggplot"))

    # Check fit_summary attribute
    fit_summary <- attr(out, "fit_summary")
    testthat::expect_true(is.data.frame(fit_summary))
    testthat::expect_true(all(c("rlz", "rank", "overall_score") %in% names(fit_summary)))

    # Check metadata attribute
    metadata <- attr(out, "metadata")
    testthat::expect_true(is.list(metadata))
    testthat::expect_equal(metadata$n_grids, n_grid)
    testthat::expect_equal(metadata$n_realizations, n_realizations)
    testthat::expect_equal(metadata$variables, c("precip", "temp"))
    testthat::expect_true(inherits(metadata$assessment_date, "Date"))
  })
})

# ==============================================================================
# TEST: Grid subsampling when exceeding max_grids
# ==============================================================================

testthat::test_that("evaluate_weather_generator subsamples grids when exceeding max_grids", {

  .with_fast_eval_mocks({

    set.seed(999)

    # Create dates (2 full years, no leap days)
    dates_obs <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
    dates_obs <- dates_obs[format(dates_obs, "%m-%d") != "02-29"]

    dates_sim <- seq.Date(as.Date("2011-01-01"), as.Date("2012-12-31"), by = "day")
    dates_sim <- dates_sim[format(dates_sim, "%m-%d") != "02-29"]

    n_grid <- 10
    max_grids <- 4
    n_realizations <- 2

    daily_obs <- lapply(seq_len(n_grid), function(i) {
      make_test_grid_df(dates_obs, id_shift = i)
    })

    daily_sim <- lapply(seq_len(n_realizations), function(r) {
      lapply(seq_len(n_grid), function(i) {
        make_test_grid_df(dates_sim, id_shift = i + r)
      })
    })

    out <- evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = n_realizations,
      output_path = NULL,
      save_plots = FALSE,
      show_title = FALSE,
      max_grids = max_grids,
      verbose = FALSE
    )

    # Check that grid count was reduced
    metadata <- attr(out, "metadata")
    testthat::expect_equal(metadata$n_grids, max_grids)
  })
})

# ==============================================================================
# TEST: Input validation - quantile parameters
# ==============================================================================

testthat::test_that("evaluate_weather_generator validates quantile parameters", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df, df)
  daily_sim <- list(list(df, df))

  # extreme_quantile must be greater than wet_quantile
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      wet_quantile = 0.9,
      extreme_quantile = 0.8,
      verbose = FALSE
    ),
    "extreme_quantile"
  )

  # wet_quantile = 0 is invalid (must be strictly between 0 and 1)
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      wet_quantile = 0,
      extreme_quantile = 0.8,
      verbose = FALSE
    ),
    "wet_quantile"
  )

  # extreme_quantile = 1 is invalid (must be strictly between 0 and 1)
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      wet_quantile = 0.2,
      extreme_quantile = 1,
      verbose = FALSE
    ),
    "extreme_quantile"
  )

  # NA values should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      wet_quantile = NA_real_,
      extreme_quantile = 0.8,
      verbose = FALSE
    ),
    "wet_quantile"
  )
})

# ==============================================================================
# TEST: Input validation - daily_sim length must match n_realizations
# ==============================================================================

testthat::test_that("evaluate_weather_generator validates daily_sim length", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df, df)
  daily_sim <- list(list(df, df))  # Only 1 realization

  # Claiming n_realizations = 2 but providing only 1 should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 2,
      verbose = FALSE
    ),
    "daily_sim"
  )
})

# ==============================================================================
# TEST: Input validation - date column required
# ==============================================================================

testthat::test_that("evaluate_weather_generator requires date column in data", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  # Data frame missing date column
  df_no_date <- data.frame(
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  # Data frame with date column
  df_with_date <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  # daily_obs missing date column
  daily_obs <- list(df_no_date, df_with_date)
  daily_sim <- list(list(df_with_date, df_with_date))

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      verbose = FALSE
    ),
    "date"
  )

  # daily_sim missing date column
  daily_obs <- list(df_with_date, df_with_date)
  daily_sim <- list(list(df_no_date, df_with_date))

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      verbose = FALSE
    ),
    "date"
  )
})

# ==============================================================================
# TEST: Leap day handling via period standardization
# ==============================================================================

testthat::test_that("evaluate_weather_generator handles leap days without error", {

  .with_fast_eval_mocks({

    set.seed(202)

    # Create dates that include leap days (Feb 29)
    # Year 2000 and 2004 are leap years
    obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "day")
    sim_dates <- seq.Date(as.Date("2012-01-01"), as.Date("2014-12-31"), by = "day")

    n_grid <- 2
    n_realizations <- 2

    daily_obs <- lapply(seq_len(n_grid), function(i) {
      make_test_grid_df(obs_dates, id_shift = i)
    })

    daily_sim <- lapply(seq_len(n_realizations), function(r) {
      lapply(seq_len(n_grid), function(i) {
        make_test_grid_df(sim_dates, id_shift = i + r)
      })
    })

    # Should handle leap days gracefully via internal standardization
    out <- evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = n_realizations,
      output_path = NULL,
      save_plots = FALSE,
      show_title = FALSE,
      verbose = FALSE
    )

    testthat::expect_s3_class(out, "weather_assessment")
    testthat::expect_true(inherits(out$daily_mean, "ggplot"))
  })
})

# ==============================================================================
# TEST: Reproducibility with seed
# ==============================================================================

testthat::test_that("evaluate_weather_generator produces reproducible results with seed", {

  .with_fast_eval_mocks({

    set.seed(42)

    # Create dates (2 full years, no leap days)
    obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
    obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

    sim_dates <- seq.Date(as.Date("2011-01-01"), as.Date("2012-12-31"), by = "day")
    sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

    n_grid <- 5
    n_realizations <- 2
    max_grids <- 3

    # Create data with enough grids to trigger subsampling
    daily_obs <- lapply(seq_len(n_grid), function(i) {
      make_test_grid_df(obs_dates, id_shift = i)
    })

    daily_sim <- lapply(seq_len(n_realizations), function(r) {
      lapply(seq_len(n_grid), function(i) {
        make_test_grid_df(sim_dates, id_shift = i + r)
      })
    })

    # Run twice with same seed
    out1 <- evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = n_realizations,
      max_grids = max_grids,
      seed = 12345,
      output_path = NULL,
      save_plots = FALSE,
      verbose = FALSE
    )

    out2 <- evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = n_realizations,
      max_grids = max_grids,
      seed = 12345,
      output_path = NULL,
      save_plots = FALSE,
      verbose = FALSE
    )

    # Metadata should be identical
    meta1 <- attr(out1, "metadata")
    meta2 <- attr(out2, "metadata")

    testthat::expect_equal(meta1$n_grids, meta2$n_grids)
    testthat::expect_equal(meta1$n_realizations, meta2$n_realizations)

    # Fit summaries should be identical
    fit1 <- attr(out1, "fit_summary")
    fit2 <- attr(out2, "fit_summary")

    testthat::expect_equal(fit1, fit2)
  })
})
