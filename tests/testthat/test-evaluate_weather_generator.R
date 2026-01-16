# ==============================================================================
# Tests for evaluate_weather_generator()
# ==============================================================================
#
# This test file validates the stochastic weather generator evaluation function,
# including:
#
#   - Output structure and class attributes
#   - Grid subsampling behavior
#   - Input validation for all parameters
#   - Leap day handling via period standardization
#
# All tests use synthetic data to ensure reproducibility and isolation.
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


# ==============================================================================
# TEST: Basic output structure and class
# ==============================================================================

testthat::test_that("evaluate_weather_generator returns weather_assessment with expected structure", {


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


# ==============================================================================
# TEST: Grid subsampling when exceeding max_grids
# ==============================================================================

testthat::test_that("evaluate_weather_generator subsamples grids when exceeding max_grids", {

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


# ==============================================================================
# TEST: Input validation - max_grids parameter
# ==============================================================================

testthat::test_that("evaluate_weather_generator rejects invalid max_grids values", {

  # Helper to create minimal valid inputs
  make_minimal_inputs <- function() {
    dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
    dates <- dates[format(dates, "%m-%d") != "02-29"]

    df <- data.frame(
      date = dates,
      precip = rep(1, length(dates)),
      temp = rep(10, length(dates))
    )

    list(
      daily_obs = list(df, df),
      daily_sim = list(list(df, df))
    )
  }

  dat <- make_minimal_inputs()

  # max_grids = 0 should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = dat$daily_sim,
      daily_obs = dat$daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      max_grids = 0,
      verbose = FALSE
  ),
    "max_grids"
  )

  # max_grids = NA should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = dat$daily_sim,
      daily_obs = dat$daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      max_grids = NA_real_,
      verbose = FALSE
  ),
    "max_grids"
  )

  # max_grids as non-integer should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = dat$daily_sim,
      daily_obs = dat$daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      max_grids = 2.5,
      verbose = FALSE
  ),
    "max_grids"
  )
})


# ==============================================================================
# TEST: Input validation - variables must include precip
# ==============================================================================

testthat::test_that("evaluate_weather_generator requires precip in variables", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df, df)
  daily_sim <- list(list(df, df))

  # variables without precip should error

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("temp"),
      n_realizations = 1,
      verbose = FALSE
  ),
    "precip"
  )
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
# TEST: Input validation - seed parameter
# ==============================================================================

testthat::test_that("evaluate_weather_generator validates seed parameter", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df)
  daily_sim <- list(list(df))

  # seed as non-integer numeric should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      seed = 1.5,
      verbose = FALSE
  ),
    "seed"
  )

  # seed as character should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      seed = "abc",
      verbose = FALSE
  ),
    "seed"
  )
})


# ==============================================================================
# TEST: Input validation - verbose parameter
# ==============================================================================

testthat::test_that("evaluate_weather_generator validates verbose parameter", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df)
  daily_sim <- list(list(df))

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      output_path = NULL,
      save_plots = FALSE,
      verbose = "yes"
  ),
    "verbose"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      output_path = NULL,
      save_plots = FALSE,
      verbose = c(TRUE, FALSE)
  ),
    "verbose"
  )
})


# ==============================================================================
# TEST: Leap day handling via period standardization
# ==============================================================================

testthat::test_that("evaluate_weather_generator handles leap days without error", {

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


# ==============================================================================
# TEST: NULL input validation
# ==============================================================================

testthat::test_that("evaluate_weather_generator rejects NULL required inputs", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df)
  daily_sim <- list(list(df))

  # daily_sim = NULL
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = NULL,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1
  ),
    "daily_sim"
  )

  # daily_obs = NULL
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = NULL,
      variables = c("precip", "temp"),
      n_realizations = 1,
      verbose = FALSE
  ),
    "daily_obs"
  )

  # variables = NULL
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = NULL,
      n_realizations = 1,
      verbose = FALSE
  ),
    "variables"
  )

  # n_realizations = NULL
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = NULL,
      verbose = FALSE
  ),
    "n_realizations"
  )
})


# ==============================================================================
# TEST: Variables not found in daily_obs
# ==============================================================================

testthat::test_that("evaluate_weather_generator errors when variables not in daily_obs", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df)
  daily_sim <- list(list(df))

  # Request variable that does not exist
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "humidity"),
      n_realizations = 1,
      verbose = FALSE
  ),
    "humidity"
  )
})


# ==============================================================================
# TEST: Empty daily_obs list
# ==============================================================================

testthat::test_that("evaluate_weather_generator rejects empty daily_obs", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list()
  daily_sim <- list(list(df))

  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1,
      verbose = FALSE
  ),
    "daily_obs"
  )
})


# ==============================================================================
# TEST: n_realizations validation
# ==============================================================================

testthat::test_that("evaluate_weather_generator validates n_realizations parameter", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(
    date = dates,
    precip = rep(1, length(dates)),
    temp = rep(10, length(dates))
  )

  daily_obs <- list(df)
  daily_sim <- list(list(df))

  # n_realizations = 0 should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 0,
      verbose = FALSE
  ),
    "n_realizations"
  )

  # n_realizations as non-integer should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = 1.5,
      verbose = FALSE
  ),
    "n_realizations"
  )

  # n_realizations as negative should error
  testthat::expect_error(
    evaluate_weather_generator(
      daily_sim = daily_sim,
      daily_obs = daily_obs,
      variables = c("precip", "temp"),
      n_realizations = -1,
      verbose = FALSE
  ),
    "n_realizations"
  )
})


# ==============================================================================
# TEST: Reproducibility with seed
# ==============================================================================

testthat::test_that("evaluate_weather_generator produces reproducible results with seed", {

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



