# tests/testthat/test-evaluate_weather_generator.R

testthat::test_that("evaluate_weather_generator returns weather_assessment with expected structure", {

  set.seed(123)

  # -----------------------------
  # Build small but valid dataset
  # -----------------------------
  make_grid_df <- function(dates, id_shift = 0) {
    n <- length(dates)
    data.frame(
      date = dates,
      precip = rgamma(n, shape = 2, scale = 2) + id_shift * 0.1,
      temp = rnorm(n, mean = 10 + id_shift, sd = 3)
    )
  }

  # 3 full years (no leap days) -> stable and fast
  obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2003-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  sim_dates <- seq.Date(as.Date("2011-01-01"), as.Date("2013-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  n_grids <- 3
  n_rlz <- 2

  daily.obs <- lapply(seq_len(n_grids), function(i) make_grid_df(obs_dates, id_shift = i))

  daily.sim <- lapply(seq_len(n_rlz), function(r) {
    lapply(seq_len(n_grids), function(i) make_grid_df(sim_dates, id_shift = i + r))
  })

  out <- evaluate_weather_generator(
    daily.sim = daily.sim,
    daily.obs = daily.obs,
    variables = c("precip", "temp"),
    realization.num = n_rlz,
    wet.quantile = 0.2,
    extreme.quantile = 0.8,
    output.path = NULL,   # force no disk writes
    save.plots = TRUE,    # will be overridden to FALSE internally
    show.title = FALSE,
    max.grids = 25
  )

  testthat::expect_s3_class(out, "weather_assessment")
  testthat::expect_true(is.list(out))
  testthat::expect_true(length(out) > 0)

  # plots should be ggplot objects (at least the canonical ones)
  testthat::expect_true(inherits(out$daily_mean, "ggplot"))
  testthat::expect_true(inherits(out$daily_sd, "ggplot"))
  testthat::expect_true(inherits(out$annual_precip, "ggplot"))

  # attributes expected by the function contract
  fit_summary <- attr(out, "fit_summary")
  metadata <- attr(out, "metadata")

  testthat::expect_true(is.data.frame(fit_summary))
  testthat::expect_true(all(c("rlz", "rank", "overall_score") %in% names(fit_summary)))
  testthat::expect_true(is.list(metadata))
  testthat::expect_equal(metadata$n_grids, n_grids)
  testthat::expect_equal(metadata$realization.num, n_rlz)
  testthat::expect_equal(metadata$variables, c("precip", "temp"))
  testthat::expect_true(inherits(metadata$assessment.date, "Date"))
})

testthat::test_that("evaluate_weather_generator subsamples grids when exceeding max.grids", {

  set.seed(999)

  make_grid_df <- function(dates, id_shift = 0) {
    n <- length(dates)
    data.frame(
      date = dates,
      precip = rgamma(n, shape = 2, scale = 1) + id_shift * 0.05,
      temp = rnorm(n, mean = 12 + id_shift, sd = 2)
    )
  }

  dates_obs <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates_obs <- dates_obs[format(dates_obs, "%m-%d") != "02-29"]

  dates_sim <- seq.Date(as.Date("2011-01-01"), as.Date("2012-12-31"), by = "day")
  dates_sim <- dates_sim[format(dates_sim, "%m-%d") != "02-29"]

  n_grids <- 10
  max_grids <- 4
  n_rlz <- 2

  daily.obs <- lapply(seq_len(n_grids), function(i) make_grid_df(dates_obs, i))
  daily.sim <- lapply(seq_len(n_rlz), function(r) {
    lapply(seq_len(n_grids), function(i) make_grid_df(dates_sim, i + r))
  })

  out <- evaluate_weather_generator(
    daily.sim = daily.sim,
    daily.obs = daily.obs,
    variables = c("precip", "temp"),
    realization.num = n_rlz,
    output.path = NULL,
    show.title = FALSE,
    max.grids = max_grids
  )

  metadata <- attr(out, "metadata")
  testthat::expect_equal(metadata$n_grids, max_grids)
})

testthat::test_that("evaluate_weather_generator rejects invalid max.grids", {

  make_minimal <- function() {
    dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
    dates <- dates[format(dates, "%m-%d") != "02-29"]
    df <- data.frame(date = dates, precip = rep(1, length(dates)), temp = rep(10, length(dates)))
    daily.obs <- list(df, df)
    daily.sim <- list(list(df, df))
    list(daily.obs = daily.obs, daily.sim = daily.sim)
  }

  dat <- make_minimal()

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = dat$daily.sim,
      daily.obs = dat$daily.obs,
      variables = c("precip", "temp"),
      realization.num = 1,
      max.grids = 0
    ),
    "max.grids"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = dat$daily.sim,
      daily.obs = dat$daily.obs,
      variables = c("precip", "temp"),
      realization.num = 1,
      max.grids = NA_real_
    ),
    "max.grids"
  )
})

testthat::test_that("evaluate_weather_generator input validation: variables must include precip", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(date = dates, precip = rep(1, length(dates)), temp = rep(10, length(dates)))
  daily.obs <- list(df, df)
  daily.sim <- list(list(df, df))

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = daily.sim,
      daily.obs = daily.obs,
      variables = c("temp"),
      realization.num = 1
    ),
    "must include 'precip'"
  )
})

testthat::test_that("evaluate_weather_generator input validation: quantiles must be valid and ordered", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(date = dates, precip = rep(1, length(dates)), temp = rep(10, length(dates)))
  daily.obs <- list(df, df)
  daily.sim <- list(list(df, df))

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = daily.sim,
      daily.obs = daily.obs,
      variables = c("precip", "temp"),
      realization.num = 1,
      wet.quantile = 0.9,
      extreme.quantile = 0.8
    ),
    "extreme.quantile.*greater"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = daily.sim,
      daily.obs = daily.obs,
      variables = c("precip", "temp"),
      realization.num = 1,
      wet.quantile = 0,
      extreme.quantile = 0.8
    ),
    "wet.quantile"
  )

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = daily.sim,
      daily.obs = daily.obs,
      variables = c("precip", "temp"),
      realization.num = 1,
      wet.quantile = 0.2,
      extreme.quantile = 1
    ),
    "extreme.quantile"
  )
})

testthat::test_that("evaluate_weather_generator input validation: daily.sim length must match realization.num", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df <- data.frame(date = dates, precip = rep(1, length(dates)), temp = rep(10, length(dates)))
  daily.obs <- list(df, df)
  daily.sim <- list(list(df, df))  # only 1 realization

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = daily.sim,
      daily.obs = daily.obs,
      variables = c("precip", "temp"),
      realization.num = 2
    ),
    "Length of 'daily.sim' must equal 'realization.num'"
  )
})

testthat::test_that("evaluate_weather_generator rejects missing date column", {

  dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  dates <- dates[format(dates, "%m-%d") != "02-29"]

  df_bad <- data.frame(precip = rep(1, length(dates)), temp = rep(10, length(dates)))
  df_ok  <- data.frame(date = dates, precip = rep(1, length(dates)), temp = rep(10, length(dates)))

  daily.obs <- list(df_bad, df_ok)
  daily.sim <- list(list(df_ok, df_ok))

  testthat::expect_error(
    evaluate_weather_generator(
      daily.sim = daily.sim,
      daily.obs = daily.obs,
      variables = c("precip", "temp"),
      realization.num = 1
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
      precip = rgamma(n, 2, 2) + id_shift * 0.05,
      temp = rnorm(n, 10 + id_shift, 3)
    )
  }

  # includes 2000-02-29
  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "day")
  sim_dates <- seq.Date(as.Date("2012-01-01"), as.Date("2014-12-31"), by = "day")

  n_grids <- 2
  n_rlz <- 2

  daily.obs <- lapply(seq_len(n_grids), function(i) make_grid_df(obs_dates, i))
  daily.sim <- lapply(seq_len(n_rlz), function(r) {
    lapply(seq_len(n_grids), function(i) make_grid_df(sim_dates, i + r))
  })

  out <- evaluate_weather_generator(
    daily.sim = daily.sim,
    daily.obs = daily.obs,
    variables = c("precip", "temp"),
    realization.num = n_rlz,
    output.path = NULL,
    show.title = FALSE
  )

  testthat::expect_s3_class(out, "weather_assessment")
  testthat::expect_true(inherits(out$daily_mean, "ggplot"))
})
