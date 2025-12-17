testthat::test_that("resample_weather_dates returns Date vector of correct length", {

  set.seed(123)

  ## ---- observed data: 6 years (safe for KNN) ----
  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  dates.d <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  PRCP <- rgamma(length(obs_dates), shape = 2, scale = 2)
  TEMP <- rnorm(length(obs_dates), mean = 10, sd = 3)

  ## ---- simulated dates: 1 year ----
  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    PRCP_FINAL_ANNUAL_SIM = 500,
    ANNUAL_PRCP = tapply(PRCP, dates.d$wyear, sum),
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2020,
    k1 = 1,
    ymax = 1,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    knn.annual.sample.num = 3,
    month.start = 1,
    seed = 123
  )

  testthat::expect_s3_class(out, "Date")
  testthat::expect_length(out, nrow(sim.dates.d))
  testthat::expect_false(anyNA(out))
})

# ------------------------------------------------------------
# Calendar-year logic: forbid cross-year observed transitions
# ------------------------------------------------------------
testthat::test_that("calendar-year mode forbids observed Dec->Jan transitions", {

  set.seed(42)

  obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  dates.d <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  PRCP <- seq_along(obs_dates)
  TEMP <- seq_along(obs_dates)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = as.integer(format(sim_dates, "%Y"))
  )

  out <- resample_weather_dates(
    PRCP_FINAL_ANNUAL_SIM = 100,
    ANNUAL_PRCP = c(100, 100),
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2020,
    k1 = 1,
    ymax = 1,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    month.start = 1,
    seed = 999
  )

  obs_year <- dates.d$wyear[match(out, dates.d$date)]
  obs_mon  <- as.integer(format(out, "%m"))
  obs_day  <- as.integer(format(out, "%d"))

  dec31 <- which(obs_mon == 12 & obs_day == 31)
  if (length(dec31) > 0 && max(dec31) < length(out)) {
    testthat::expect_true(all(obs_year[dec31 + 1] == obs_year[dec31]))
  }
})

testthat::test_that("water-year mode allows observed Dec->Jan transitions", {

  set.seed(7)

  ## ---- observed data: 6 water years ----
  obs_dates <- seq.Date(as.Date("2000-10-01"), as.Date("2006-09-30"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  dates.d <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = ifelse(
      as.integer(format(obs_dates, "%m")) >= 10,
      as.integer(format(obs_dates, "%Y")) + 1,
      as.integer(format(obs_dates, "%Y"))
    )
  )

  PRCP <- rgamma(length(obs_dates), 2, 2)
  TEMP <- rnorm(length(obs_dates), 10, 3)

  sim_dates <- seq.Date(as.Date("2020-10-01"), as.Date("2021-09-30"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2021
  )

  out <- resample_weather_dates(
    PRCP_FINAL_ANNUAL_SIM = 600,
    ANNUAL_PRCP = tapply(PRCP, dates.d$wyear, sum),
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2020,
    k1 = 1,
    ymax = 1,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    knn.annual.sample.num = 3,
    month.start = 10,
    seed = 321
  )

  mons <- as.integer(format(out, "%m"))
  testthat::expect_true(any(mons[-1] == 1 & mons[-length(mons)] == 12))
})

# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------
testthat::test_that("resample_weather_dates is reproducible with same seed", {

  set.seed(123)

  n <- 5 * 365
  dates_obs <- seq.Date(as.Date("2000-01-01"), by = "day", length.out = n)

  dates.d <- data.frame(
    date  = dates_obs,
    month = as.integer(format(dates_obs, "%m")),
    day   = as.integer(format(dates_obs, "%d")),
    wyear = as.integer(format(dates_obs, "%Y"))
  )

  PRCP <- rgamma(n, 2, 2)
  TEMP <- rnorm(n, 15, 5)
  ANNUAL_PRCP <- tapply(PRCP, dates.d$wyear, sum)

  sim_dates <- seq.Date(as.Date("2010-01-01"), by = "day", length.out = 365)
  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2010
  )

  out1 <- resample_weather_dates(
    mean(ANNUAL_PRCP),
    ANNUAL_PRCP,
    PRCP,
    TEMP,
    START_YEAR_SIM = 2010,
    k1 = 1,
    ymax = 1,
    dates.d,
    sim.dates.d,
    seed = 42
  )

  out2 <- resample_weather_dates(
    mean(ANNUAL_PRCP),
    ANNUAL_PRCP,
    PRCP,
    TEMP,
    START_YEAR_SIM = 2010,
    k1 = 1,
    ymax = 1,
    dates.d,
    sim.dates.d,
    seed = 42
  )

  testthat::expect_identical(out1, out2)
})

# ------------------------------------------------------------
# Input validation
# ------------------------------------------------------------
testthat::test_that("resample_weather_dates rejects invalid month.start", {

  expect_error(
    resample_weather_dates(
      PRCP_FINAL_ANNUAL_SIM = 100,
      ANNUAL_PRCP = 100,
      PRCP = 1,
      TEMP = 1,
      START_YEAR_SIM = 2000,
      k1 = 1,
      ymax = 1,
      dates.d = data.frame(date = Sys.Date(), month = 1, day = 1, wyear = 2000),
      sim.dates.d = data.frame(month = 1, day = 1, wyear = 2000),
      month.start = 13
    ),
    "month.start must be in 1:12"
  )
})

# ------------------------------------------------------------
# Markov integration
# ------------------------------------------------------------
testthat::test_that("resample_weather_dates produces non-degenerate sequence", {

  set.seed(99)

  dates_obs <- seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3 * 365)
  dates.d <- data.frame(
    date  = dates_obs,
    month = as.integer(format(dates_obs, "%m")),
    day   = as.integer(format(dates_obs, "%d")),
    wyear = as.integer(format(dates_obs, "%Y"))
  )

  PRCP <- rgamma(length(dates_obs), 2, 2)
  TEMP <- rnorm(length(dates_obs))
  ANNUAL_PRCP <- tapply(PRCP, dates.d$wyear, sum)

  sim_dates <- seq.Date(as.Date("2010-01-01"), by = "day", length.out = 365)
  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2010
  )

  out <- resample_weather_dates(
    mean(ANNUAL_PRCP),
    ANNUAL_PRCP,
    PRCP,
    TEMP,
    START_YEAR_SIM = 2010,
    k1 = 1,
    ymax = 1,
    dates.d,
    sim.dates.d,
    seed = 1
  )

  testthat::expect_gt(length(unique(out)), 30)
})

testthat::test_that("resample_weather_dates does not fail with small candidate windows", {

  set.seed(101)

  obs_dates <- seq.Date(as.Date("2000-01-01"), by = "day", length.out = 6 * 365)

  dates.d <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  PRCP <- rgamma(length(obs_dates), 2, 1)
  TEMP <- rnorm(length(obs_dates), 10, 2)
  ANNUAL_PRCP <- tapply(PRCP, dates.d$wyear, sum)

  sim_dates <- seq.Date(as.Date("2010-01-01"), by = "day", length.out = 365)
  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2010
  )

  expect_error(
    resample_weather_dates(
      PRCP_FINAL_ANNUAL_SIM = mean(ANNUAL_PRCP),
      ANNUAL_PRCP = ANNUAL_PRCP,
      PRCP = PRCP,
      TEMP = TEMP,
      START_YEAR_SIM = 2010,
      k1 = 1,
      ymax = 1,
      dates.d = dates.d,
      sim.dates.d = sim.dates.d,
      knn.annual.sample.num = 2,
      seed = 1
    ),
    NA
  )
})
