
test_that("resample_weather_dates returns Date vector of correct length", {

  set.seed(123)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- rgamma(length(obs_dates), shape = 2, scale = 2)
  obs_daily_temp   <- rnorm(length(obs_dates), mean = 10, sd = 3)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    sim_annual_precip  = 500,
    obs_annual_precip  = tapply(obs_daily_precip, obs_dates_df$wyear, sum),
    obs_daily_precip   = obs_daily_precip,
    obs_daily_temp     = obs_daily_temp,
    sim_year_start     = 2020,
    k1                 = 1,
    n_sim_years        = 1,
    obs_dates_df       = obs_dates_df,
    sim_dates_df       = sim_dates_df,
    knn.annual.sample.num = 3,
    month.start        = 1,
    seed               = 123
  )

  testthat::expect_s3_class(out, "Date")
  testthat::expect_length(out, nrow(sim_dates_df))
  testthat::expect_false(anyNA(out))
})

testthat::test_that("calendar-year mode forbids observed Dec->Jan transitions", {

  set.seed(42)

  obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- seq_along(obs_dates)
  obs_daily_temp   <- seq_along(obs_dates)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    sim_annual_precip = 100,
    obs_annual_precip = c(100, 100),
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    sim_year_start    = 2020,
    k1                = 1,
    n_sim_years       = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    month.start       = 1,
    seed              = 999
  )

  obs_year <- obs_dates_df$wyear[match(out, obs_dates_df$date)]
  obs_mon  <- as.integer(format(out, "%m"))
  obs_day  <- as.integer(format(out, "%d"))

  dec31 <- which(obs_mon == 12 & obs_day == 31)
  if (length(dec31) > 0 && max(dec31) < length(out)) {
    testthat::expect_true(all(obs_year[dec31 + 1] == obs_year[dec31]))
  }
})

testthat::test_that("water-year mode allows observed Dec->Jan transitions", {

  set.seed(7)

  obs_dates <- seq.Date(as.Date("2000-10-01"), as.Date("2006-09-30"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = ifelse(
      as.integer(format(obs_dates, "%m")) >= 10,
      as.integer(format(obs_dates, "%Y")) + 1,
      as.integer(format(obs_dates, "%Y"))
    )
  )

  obs_daily_precip <- rgamma(length(obs_dates), 2, 2)
  obs_daily_temp   <- rnorm(length(obs_dates), 10, 3)

  sim_dates <- seq.Date(as.Date("2020-10-01"), as.Date("2021-09-30"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2021
  )

  out <- resample_weather_dates(
    sim_annual_precip = 600,
    obs_annual_precip = tapply(obs_daily_precip, obs_dates_df$wyear, sum),
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    sim_year_start    = 2020,
    k1                = 1,
    n_sim_years       = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    knn.annual.sample.num = 3,
    month.start       = 10,
    seed              = 321
  )

  mons <- as.integer(format(out, "%m"))
  testthat::expect_true(any(mons[-1] == 1 & mons[-length(mons)] == 12))
})

testthat::test_that("resample_weather_dates is reproducible with same seed", {

  n <- 5 * 365
  dates_obs <- seq.Date(as.Date("2000-01-01"), by = "day", length.out = n)

  obs_dates_df <- data.frame(
    date  = dates_obs,
    month = as.integer(format(dates_obs, "%m")),
    day   = as.integer(format(dates_obs, "%d")),
    wyear = as.integer(format(dates_obs, "%Y"))
  )

  obs_daily_precip <- rgamma(n, 2, 2)
  obs_daily_temp   <- rnorm(n, 15, 5)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2010-01-01"), by = "day", length.out = 365)
  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2010
  )

  out1 <- resample_weather_dates(
    mean(obs_annual_precip),
    obs_annual_precip,
    obs_daily_precip,
    obs_daily_temp,
    sim_year_start = 2010,
    k1             = 1,
    n_sim_years    = 1,
    obs_dates_df,
    sim_dates_df,
    seed = 42
  )

  out2 <- resample_weather_dates(
    mean(obs_annual_precip),
    obs_annual_precip,
    obs_daily_precip,
    obs_daily_temp,
    sim_year_start = 2010,
    k1             = 1,
    n_sim_years    = 1,
    obs_dates_df,
    sim_dates_df,
    seed = 42
  )

  testthat::expect_identical(out1, out2)
})

testthat::test_that("resample_weather_dates rejects invalid month.start", {

  expect_error(
    resample_weather_dates(
      sim_annual_precip = 100,
      obs_annual_precip = 100,
      obs_daily_precip  = 1,
      obs_daily_temp    = 1,
      sim_year_start    = 2000,
      k1                = 1,
      n_sim_years       = 1,
      obs_dates_df      = data.frame(date = Sys.Date(), month = 1, day = 1, wyear = 2000),
      sim_dates_df      = data.frame(month = 1, day = 1, wyear = 2000),
      month.start       = 13
    ),
    "month.start must be a single integer between 1 and 12"
  )
})




