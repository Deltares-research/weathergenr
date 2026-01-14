testthat::test_that("evaluate_weather_generator returns weather_assessment with expected structure", {

  set.seed(123)

  make_grid_df <- function(dates, id_shift = 0) {
    n <- length(dates)
    data.frame(
      date = dates,
      prcp = rgamma(n, shape = 2, scale = 2) + id_shift * 0.1,
      temp = rnorm(n, mean = 10 + id_shift, sd = 3)
    )
  }

  obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2003-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  sim_dates <- seq.Date(as.Date("2011-01-01"), as.Date("2013-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  n_grids <- 3
  n_realizations <- 2

  obs_data <- lapply(seq_len(n_grids), function(i) make_grid_df(obs_dates, id_shift = i))

  sim_data <- lapply(seq_len(n_realizations), function(r) {
    lapply(seq_len(n_grids), function(i) make_grid_df(sim_dates, id_shift = i + r))
  })

  out <- evaluate_weather_generator(
    sim_data = sim_data,
    obs_data = obs_data,
    vars = c("prcp", "temp"),
    n_realizations = n_realizations,
    wet_q = 0.2,
    extreme_q = 0.8,
    out_dir = NULL,     # no disk writes
    save_plot = TRUE,  # function should force FALSE if out_dir is NULL
    show_title = FALSE,
    max_grid = 25
  )

  testthat::expect_s3_class(out, "weather_assessment")
  testthat::expect_true(is.list(out))
  testthat::expect_true(length(out) > 0)

  testthat::expect_true(inherits(out$daily_mean, "ggplot"))
  testthat::expect_true(inherits(out$daily_sd, "ggplot"))
  testthat::expect_true(inherits(out$annual_prcp, "ggplot"))

  fit_summary <- attr(out, "fit_summary")
  metadata <- attr(out, "metadata")

  testthat::expect_true(is.data.frame(fit_summary))
  testthat::expect_true(all(c("rlz", "rank", "overall_score") %in% names(fit_summary)))

  testthat::expect_true(is.list(metadata))
  testthat::expect_equal(metadata$n_grids, n_grids)
  testthat::expect_equal(metadata$n_realizations, n_realizations)
  testthat::expect_equal(metadata$vars, c("prcp", "temp"))
  testthat::expect_true(inherits(metadata$assessment_date, "Date"))
})

testthat::test_that("evaluate_weather_generator subsamples grids when exceeding max_grid", {

  set.seed(999)

  make_grid_df <- function(dates, id_shift = 0) {
    n <- length(dates)
    data.frame(
      date = dates,
      prcp = rgamma(n, shape = 2, scale = 1) + id_shift * 0.05,
      temp = rnorm(n, mean = 12 + id_shift, sd = 2)
    )
  }

  dates_obs <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates_obs <- dates_obs[format(dates_obs, "%m-%d") != "02-29"]

  dates_sim <- seq.Date(as.Date("2011-01-01"), as.Date("2012-12-31"), by = "day")
  dates_sim <- dates_sim[format(dates_sim, "%m-%d") != "02-29"]

  n_grids <- 10
  max_grid <- 4
  n_realizations <- 2

  obs_data <- lapply(seq_len(n_grids), function(i) make_grid_df(dates_obs, i))
  sim_data <- lapply(seq_len(n_realizations), function(r) {
    lapply(seq_len(n_grids), function(i) make_grid_df(dates_sim, i + r))
  })

  out <- evaluate_weather_generator(
    sim_data = sim_data,
    obs_data = obs_data,
    vars = c("prcp", "temp"),
    n_realizations = n_realizations,
    out_dir = NULL,
    show_title = FALSE,
    max_grid = max_grid
  )

  metadata <- attr(out, "metadata")
  testthat::expect_equal(metadata$n_grids, max_grid)
})

testthat::test_that("evaluate_weather_generator rejects invalid max_grid", {

  make_minimal <- function() {
    dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
    dates <- dates[format(dates, "%m-%d") != "02-29"]
    df <- data.frame(date = dates, prcp = rep(1, length(dates)), temp = rep(10, length(dates)))
    obs_data <- list(df, df)
    sim_data <- list(list(df, df))
    list(obs_data = obs_data, sim_data = sim_data)
  }

  dat <- make_minimal()

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = dat$sim_data,
      obs_data = dat$obs_data,
      vars = c("prcp", "temp"),
      n_realizations = 1,
      max_grid = 0
    ),
    "max_grid"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = dat$sim_data,
      obs_data = dat$obs_data,
      vars = c("prcp", "temp"),
      n_realizations = 1,
      max_grid = NA_real_
    ),
    "max_grid"
  )
})

testthat::test_that("evaluate_weather_generator input validation: vars must include prcp", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(date = dates, prcp = rep(1, length(dates)), temp = rep(10, length(dates)))
  obs_data <- list(df, df)
  sim_data <- list(list(df, df))

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = sim_data,
      obs_data = obs_data,
      vars = c("temp"),
      n_realizations = 1
    ),
    "include 'prcp'"
  )
})

testthat::test_that("evaluate_weather_generator input validation: quantiles must be valid and ordered", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(date = dates, prcp = rep(1, length(dates)), temp = rep(10, length(dates)))
  obs_data <- list(df, df)
  sim_data <- list(list(df, df))

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = sim_data,
      obs_data = obs_data,
      vars = c("prcp", "temp"),
      n_realizations = 1,
      wet_q = 0.9,
      extreme_q = 0.8
    ),
    "extreme_q.*greater|extreme"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = sim_data,
      obs_data = obs_data,
      vars = c("prcp", "temp"),
      n_realizations = 1,
      wet_q = 0,
      extreme_q = 0.8
    ),
    "wet_q"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = sim_data,
      obs_data = obs_data,
      vars = c("prcp", "temp"),
      n_realizations = 1,
      wet_q = 0.2,
      extreme_q = 1
    ),
    "extreme_q"
  )
})

testthat::test_that("evaluate_weather_generator input validation: sim_data length must match n_realizations", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(date = dates, prcp = rep(1, length(dates)), temp = rep(10, length(dates)))
  obs_data <- list(df, df)
  sim_data <- list(list(df, df))  # 1 realization

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = sim_data,
      obs_data = obs_data,
      vars = c("prcp", "temp"),
      n_realizations = 2
    ),
    "Length.*sim_data.*n_realizations|length.*sim_data"
  )
})

testthat::test_that("evaluate_weather_generator rejects missing date column", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df_bad <- data.frame(prcp = rep(1, length(dates)), temp = rep(10, length(dates)))
  df_ok  <- data.frame(date = dates, prcp = rep(1, length(dates)), temp = rep(10, length(dates)))

  obs_data <- list(df_bad, df_ok)
  sim_data <- list(list(df_ok, df_ok))

  testthat::expect_error(
    evaluate_weather_generator(
      sim_data = sim_data,
      obs_data = obs_data,
      vars = c("prcp", "temp"),
      n_realizations = 1
    ),
    "date"
  )
})

testthat::test_that("evaluate_weather_generator handles leap days via standardization (no error)", {

  set.seed(202)

  make_grid_df <- function(dates, id_shift = 0) {
    n <- length(dates)
    data.frame(
      date = dates,
      prcp = rgamma(n, 2, 2) + id_shift * 0.05,
      temp = rnorm(n, 10 + id_shift, 3)
    )
  }

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "day")
  sim_dates <- seq.Date(as.Date("2012-01-01"), as.Date("2014-12-31"), by = "day")

  n_grids <- 2
  n_realizations <- 2

  obs_data <- lapply(seq_len(n_grids), function(i) make_grid_df(obs_dates, i))
  sim_data <- lapply(seq_len(n_realizations), function(r) {
    lapply(seq_len(n_grids), function(i) make_grid_df(sim_dates, i + r))
  })

  out <- evaluate_weather_generator(
    sim_data = sim_data,
    obs_data = obs_data,
    vars = c("prcp", "temp"),
    n_realizations = n_realizations,
    out_dir = NULL,
    show_title = FALSE
  )

  testthat::expect_s3_class(out, "weather_assessment")
  testthat::expect_true(inherits(out$daily_mean, "ggplot"))
})
