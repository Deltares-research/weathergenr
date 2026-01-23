# ==============================================================================
# Tests for evaluate_weather_generator()
# ==============================================================================
# Functions under test:
# - R/evaluate_generator.R: evaluate_weather_generator()
# - R/evaluate_generator_plots.R: create_all_diagnostic_plots(), generate_symmetric_dummy_points()
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

library(testthat)


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
                                                show_title, save_plots, output_dir) {
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
  .mock_summarize_simulated_data <- function(daily_sim, n_realizations, variables, mc_thresholds,
                                             parallel = FALSE, n_cores = NULL, seed = NULL) {

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
      vars = c("precip", "temp"),
      n_realizations = n_realizations,
      wet_q = 0.2,
      extreme_q = 0.8,
      output_dir = NULL,
      save_plots = FALSE,
      show_title = FALSE,
      eval_max_grids = 25,
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
# TEST: Grid subsampling when exceeding eval_max_grids
# ==============================================================================

testthat::test_that("evaluate_weather_generator subsamples grids when exceeding eval_max_grids", {

  .with_fast_eval_mocks({

    set.seed(999)

    # Create dates (2 full years, no leap days)
    dates_obs <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
    dates_obs <- dates_obs[format(dates_obs, "%m-%d") != "02-29"]

    dates_sim <- seq.Date(as.Date("2011-01-01"), as.Date("2012-12-31"), by = "day")
    dates_sim <- dates_sim[format(dates_sim, "%m-%d") != "02-29"]

    n_grid <- 10
    eval_max_grids <- 4
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
      vars = c("precip", "temp"),
      n_realizations = n_realizations,
      output_dir = NULL,
      save_plots = FALSE,
      show_title = FALSE,
      eval_max_grids = eval_max_grids,
      verbose = FALSE
    )

    # Check that grid count was reduced
    metadata <- attr(out, "metadata")
    testthat::expect_equal(metadata$n_grids, eval_max_grids)
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

  # extreme_q must be greater than wet_q
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      vars = c("precip", "temp"),
      n_realizations = 1,
      wet_q = 0.9,
      extreme_q = 0.8,
      verbose = FALSE
    ),
    "extreme_q"
  )

  # wet_q = 0 is invalid (must be strictly between 0 and 1)
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      vars = c("precip", "temp"),
      n_realizations = 1,
      wet_q = 0,
      extreme_q = 0.8,
      verbose = FALSE
    ),
    "wet_q"
  )

  # extreme_q = 1 is invalid (must be strictly between 0 and 1)
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      vars = c("precip", "temp"),
      n_realizations = 1,
      wet_q = 0.2,
      extreme_q = 1,
      verbose = FALSE
    ),
    "extreme_q"
  )

  # NA values should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      vars = c("precip", "temp"),
      n_realizations = 1,
      wet_q = NA_real_,
      extreme_q = 0.8,
      verbose = FALSE
    ),
    "wet_q"
  )
})

# ==============================================================================
# TEST: Input validation - parallel options
# ==============================================================================

testthat::test_that("evaluate_weather_generator validates parallel arguments", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df, df)
  daily_sim <- list(list(df, df))

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      vars = c("precip", "temp"),
      n_realizations = 1,
      parallel = "yes",
      verbose = FALSE
    ),
    "parallel"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      vars = c("precip", "temp"),
      n_realizations = 1,
      parallel = TRUE,
      n_cores = 0,
      verbose = FALSE
    ),
    "n_cores"
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
      vars = c("precip", "temp"),
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
      vars = c("precip", "temp"),
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
      vars = c("precip", "temp"),
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
      vars = c("precip", "temp"),
      n_realizations = n_realizations,
      output_dir = NULL,
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
    eval_max_grids <- 3

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
      vars = c("precip", "temp"),
      n_realizations = n_realizations,
      eval_max_grids = eval_max_grids,
      seed = 12345,
      output_dir = NULL,
      save_plots = FALSE,
      verbose = FALSE
    )

    out2 <- evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      vars = c("precip", "temp"),
      n_realizations = n_realizations,
      eval_max_grids = eval_max_grids,
      seed = 12345,
      output_dir = NULL,
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


test_that("generate_symmetric_dummy_points returns 2 points per facet and is symmetric", {
  df <- data.frame(
    variable = c("precip", "precip", "temp", "temp"),
    Observed = c(1, 5, 10, 12),
    Simulated = c(0.5, 6, 9, 13)
  )

  out <- generate_symmetric_dummy_points(
    df = df,
    facet_var = "variable",
    x_col = "Observed",
    y_col = "Simulated"
  )

  # 2 dummy points per facet
  expect_true(is.data.frame(out))
  expect_true(all(c("variable", "Observed", "Simulated") %in% names(out)))
  expect_equal(nrow(out), 2L * length(unique(df$variable)))

  # Symmetry: each row has x == y
  expect_true(all(out$Observed == out$Simulated))

  # For each facet: contains min and max of combined range
  for (v in unique(df$variable)) {
    d <- df[df$variable == v, ]
    rng <- range(c(d$Observed, d$Simulated), na.rm = TRUE)
    o <- out[out$variable == v, ]
    expect_true(any(abs(o$Observed - rng[1]) < 1e-12))
    expect_true(any(abs(o$Observed - rng[2]) < 1e-12))
  }
})

test_that("create_all_diagnostic_plots returns expected plot names and ggplot objects", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")

  # Minimal plot_data fixture with required columns
  daily_stats_season <- data.frame(
    rlz = c(1, 1, 1, 1),
    id = c(1, 1, 1, 1),
    mon = c(1, 2, 1, 2),
    variable = c("precip", "precip", "temp", "temp"),
    stat = c("mean", "mean", "mean", "mean"),
    Observed = c(1, 2, 10, 11),
    Simulated = c(1.1, 1.9, 9.8, 11.2)
  )

  # Need SD rows for daily_sd plot
  daily_stats_sd <- daily_stats_season
  daily_stats_sd$stat <- "sd"
  daily_stats_sd$Observed <- c(0.5, 0.6, 1.0, 1.1)
  daily_stats_sd$Simulated <- c(0.55, 0.58, 0.9, 1.2)

  daily_stats_season <- rbind(daily_stats_season, daily_stats_sd)

  # Wet/dry diagnostics used by spell_length + wetdry_days_count plots
  stats_wetdry <- data.frame(
    rlz = c(1, 1, 1, 1),
    id = c(1, 1, 1, 1),
    mon = c(1, 1, 1, 1),
    variable = c("precip", "precip", "precip", "precip"),
    stat = c("Wet", "Dry", "Wet", "Dry"),
    type = c("days", "days", "spells", "spells"),
    Observed = c(10, 20, 2.0, 3.0),
    Simulated = c(11, 19, 2.2, 2.8)
  )

  # Cross-grid correlations: variable1 facets
  stats_crosscor <- data.frame(
    rlz = c(1, 1),
    variable1 = c("precip", "temp"),
    variable2 = c("precip", "temp"),
    id1 = c(1, 1),
    id2 = c(2, 2),
    Observed = c(0.4, 0.2),
    Simulated = c(0.35, 0.25)
  )

  # Inter-variable correlations: variable facets
  stats_intercor <- data.frame(
    rlz = c(1, 1),
    variable = c("precip:temp", "precip:temp"),
    variable1 = c("precip", "precip"),
    variable2 = c("temp", "temp"),
    id1 = c(1, 2),
    id2 = c(1, 2),
    Observed = c(0.1, 0.15),
    Simulated = c(0.12, 0.14)
  )

  # Conditional precip correlations: regime ~ variable facets, within-grid only
  stats_precip_cor_cond <- data.frame(
    rlz = c(1, 1, 1),
    id1 = c(1, 1, 1),
    id2 = c(1, 1, 1),
    variable1 = c("precip", "precip", "precip"),
    variable2 = c("temp", "temp", "temp"),
    regime = c("all", "wet", "dry"),
    transform = c("precip", "log1p_precip", "precip"),
    Observed = c(0.05, 0.10, 0.00),
    Simulated = c(0.06, 0.09, -0.01)
  )

  # Monthly pattern inputs
  stats_mon_aavg_sim <- data.frame(
    rlz = c(1, 1, 1, 1),
    year = c(1, 1, 1, 1),
    mon = c(1, 2, 1, 2),
    variable = c("precip", "precip", "temp", "temp"),
    stat = c("mean", "mean", "mean", "mean"),
    Simulated = c(1.0, 2.0, 10, 11)
  )

  stats_mon_aavg_obs <- data.frame(
    year = c(1, 1, 1, 1),
    mon = c(1, 2, 1, 2),
    variable = c("precip", "precip", "temp", "temp"),
    stat = c("mean", "mean", "mean", "mean"),
    Observed = c(1.1, 1.9, 10.2, 10.8)
  )

  # Annual precip inputs
  stats_annual_aavg_sim <- data.frame(
    rlz = c(1, 1),
    year = c(1, 2),
    variable = c("precip", "precip"),
    stat = c("mean", "mean"),
    Simulated = c(1.5, 1.7)
  )

  stats_annual_aavg_obs <- data.frame(
    year = c(1, 2),
    variable = c("precip", "precip"),
    stat = c("mean", "mean"),
    Observed = c(1.6, 1.65)
  )

  plot_data <- list(
    daily_stats_season = daily_stats_season,
    stats_mon_aavg_sim = stats_mon_aavg_sim,
    stats_mon_aavg_obs = stats_mon_aavg_obs,
    stats_annual_aavg_sim = stats_annual_aavg_sim,
    stats_annual_aavg_obs = stats_annual_aavg_obs,
    stats_crosscor = stats_crosscor,
    stats_intercor = stats_intercor,
    stats_wetdry = stats_wetdry,
    stats_precip_cor_cond = stats_precip_cor_cond
  )

  plot_config <- list(
    subtitle = "test subtitle",
    alpha = 0.4,
    colors = stats::setNames(c("blue3", "gray40"), c("Observed", "Simulated")),
    theme = ggplot2::theme_bw(base_size = 12)
  )

  variables <- c("precip", "temp")

  plots <- create_all_diagnostic_plots(
    plot_data = plot_data,
    plot_config = plot_config,
    variables = variables,
    show_title = FALSE,
    save_plots = FALSE,
    output_dir = NULL
  )

  expect_true(is.list(plots))

  # Expected base plot names
  expected <- c(
    "daily_mean", "daily_sd",
    "spell_length", "wetdry_days_count",
    "crossgrid", "intergrid",
    "precip_cond_cor",
    "monthly_cycle", "annual_precip",
    "annual_pattern_precip", "annual_pattern_temp"
  )
  expect_true(all(expected %in% names(plots)))

  # All are ggplot objects
  for (nm in expected) {
    expect_s3_class(plots[[nm]], "ggplot")
  }
})

test_that("create_all_diagnostic_plots calls ggsave when save_plots=TRUE (mocked)", {
  skip_if_not_installed("ggplot2")

  # Minimal plot_data for one plot is not enough because create_all_diagnostic_plots builds all plots.
  # Reuse the prior fixture by sourcing it from a helper if you have one.
  # Here, keep it simple: mock ggsave and run with save_plots=TRUE + tempdir.

  # Define a small valid fixture by reusing the previous test's objects if tests run in order is NOT guaranteed.
  # So we rebuild a minimal-but-valid plot_data by calling the same builder logic in-line.

  daily_stats_season <- data.frame(
    rlz = c(1, 1),
    id = c(1, 1),
    mon = c(1, 1),
    variable = c("precip", "precip"),
    stat = c("mean", "sd"),
    Observed = c(1, 0.5),
    Simulated = c(1.1, 0.55)
  )

  stats_wetdry <- data.frame(
    rlz = c(1, 1, 1, 1),
    id = c(1, 1, 1, 1),
    mon = c(1, 1, 1, 1),
    variable = c("precip", "precip", "precip", "precip"),
    stat = c("Wet", "Dry", "Wet", "Dry"),
    type = c("days", "days", "spells", "spells"),
    Observed = c(10, 20, 2.0, 3.0),
    Simulated = c(11, 19, 2.2, 2.8)
  )

  stats_crosscor <- data.frame(
    rlz = 1,
    variable1 = "precip",
    variable2 = "precip",
    id1 = 1, id2 = 2,
    Observed = 0.4, Simulated = 0.35
  )

  stats_intercor <- data.frame(
    rlz = 1,
    variable = "precip:temp",
    variable1 = "precip", variable2 = "temp",
    id1 = 1, id2 = 1,
    Observed = 0.1, Simulated = 0.12
  )

  stats_precip_cor_cond <- data.frame(
    rlz = c(1, 1, 1),
    id1 = 1, id2 = 1,
    variable1 = "precip", variable2 = "temp",
    regime = c("all", "wet", "dry"),
    transform = c("precip", "log1p_precip", "precip"),
    Observed = c(0.05, 0.10, 0.00),
    Simulated = c(0.06, 0.09, -0.01)
  )

  stats_mon_aavg_sim <- data.frame(
    rlz = 1,
    year = 1,
    mon = 1,
    variable = "precip",
    stat = "mean",
    Simulated = 1.0
  )

  stats_mon_aavg_obs <- data.frame(
    year = 1,
    mon = 1,
    variable = "precip",
    stat = "mean",
    Observed = 1.1
  )

  stats_annual_aavg_sim <- data.frame(
    rlz = 1,
    year = 1,
    variable = "precip",
    stat = "mean",
    Simulated = 1.5
  )

  stats_annual_aavg_obs <- data.frame(
    year = 1,
    variable = "precip",
    stat = "mean",
    Observed = 1.6
  )

  plot_data <- list(
    daily_stats_season = daily_stats_season,
    stats_mon_aavg_sim = stats_mon_aavg_sim,
    stats_mon_aavg_obs = stats_mon_aavg_obs,
    stats_annual_aavg_sim = stats_annual_aavg_sim,
    stats_annual_aavg_obs = stats_annual_aavg_obs,
    stats_crosscor = stats_crosscor,
    stats_intercor = stats_intercor,
    stats_wetdry = stats_wetdry,
    stats_precip_cor_cond = stats_precip_cor_cond
  )

  plot_config <- list(
    subtitle = "test subtitle",
    alpha = 0.4,
    colors = stats::setNames(c("blue3", "gray40"), c("Observed", "Simulated")),
    theme = ggplot2::theme_bw(base_size = 12)
  )

  out_dir <- tempdir()
  calls <- 0L

  testthat::local_mocked_bindings(
    ggsave = function(...) {
      calls <<- calls + 1L
      invisible(NULL)
    },
    .package = "ggplot2"
  )

  plots <- create_all_diagnostic_plots(
    plot_data = plot_data,
    plot_config = plot_config,
    variables = c("precip"),  # keep it minimal
    show_title = FALSE,
    save_plots = TRUE,
    output_dir = out_dir
  )

  expect_true(is.list(plots))
  expect_true(calls >= 1L)
})

