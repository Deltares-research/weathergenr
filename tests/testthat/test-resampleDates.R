test_that("resampleDates returns a Date vector of correct length", {

  set.seed(123)

  # ---- minimal synthetic observed data (2 years, no leap) ----
  obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  dates.d <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  PRCP <- rep(5, length(obs_dates))
  TEMP <- rep(10, length(obs_dates))

  # ---- simulated dates: 1 year ----
  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = as.integer(format(sim_dates, "%Y"))
  )

  out <- resampleDates(
    PRCP_FINAL_ANNUAL_SIM = c(100),
    ANNUAL_PRCP = c(100, 100),
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2020,
    k1 = 1,
    ymax = 1,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    month.start = 1,
    seed = 123
  )

  expect_s3_class(out, "Date")
  expect_length(out, nrow(sim.dates.d))
  expect_false(any(is.na(out)))
})

# ------------------------------------------------------------
# Calendar-year logic: forbid Dec 31 -> Jan 01 transitions
# ------------------------------------------------------------
test_that("calendar-year mode forbids cross-year transitions", {

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

  out <- resampleDates(
    PRCP_FINAL_ANNUAL_SIM = c(100),
    ANNUAL_PRCP = c(100, 100),
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2020,
    k1 = 1,
    ymax = 1,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    month.start = 1,   # calendar-year
    seed = 999
  )

  out_year <- as.integer(format(out, "%Y"))
  out_mon  <- as.integer(format(out, "%m"))
  out_day  <- as.integer(format(out, "%d"))

  # Identify Dec 31 positions
  dec31 <- which(out_mon == 12 & out_day == 31)

  if (length(dec31) > 0 && max(dec31) < length(out)) {
    # Next day after Dec 31 must stay in same year
    expect_true(all(out_year[dec31 + 1] == out_year[dec31]))
  }
})

# ------------------------------------------------------------
# Water-year logic: allow cross-calendar-year transitions
# ------------------------------------------------------------
test_that("water-year mode allows Dec to Jan transitions", {

  set.seed(7)

  obs_dates <- seq.Date(as.Date("2001-10-01"), as.Date("2003-09-30"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  dates.d <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = ifelse(format(obs_dates, "%m") %in% c("10","11","12"),
                   as.integer(format(obs_dates, "%Y")) + 1,
                   as.integer(format(obs_dates, "%Y")))
  )

  PRCP <- rep(5, length(obs_dates))
  TEMP <- rep(10, length(obs_dates))

  sim_dates <- seq.Date(as.Date("2020-10-01"), as.Date("2021-09-30"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = ifelse(format(sim_dates, "%m") %in% c("10","11","12"),
                   2021, 2021)
  )

  out <- resampleDates(
    PRCP_FINAL_ANNUAL_SIM = c(100),
    ANNUAL_PRCP = c(100, 100),
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2020,
    k1 = 1,
    ymax = 1,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    month.start = 10,   # water-year
    seed = 321
  )

  expect_s3_class(out, "Date")
  expect_false(any(is.na(out)))

  # At least one Dec->Jan transition should exist
  mons <- as.integer(format(out, "%m"))
  expect_true(any(mons[-1] == 1 & mons[-length(mons)] == 12))
})

# ------------------------------------------------------------
# First-day logic: matches simulated calendar
# ------------------------------------------------------------
test_that("first simulated day matches simulated month/day", {

  set.seed(11)

  obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  dates.d <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = 2001
  )

  PRCP <- rep(1, length(obs_dates))
  TEMP <- rep(1, length(obs_dates))

  sim_dates <- seq.Date(as.Date("2020-03-01"), as.Date("2021-02-28"), by = "day")

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resampleDates(
    PRCP_FINAL_ANNUAL_SIM = c(50),
    ANNUAL_PRCP = c(50),
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2020,
    k1 = 1,
    ymax = 1,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    month.start = 3,
    seed = 55
  )

  expect_equal(
    format(out[1], "%m-%d"),
    format(sim_dates[1], "%m-%d")
  )
})

test_that("resampleDates runs end-to-end and is reproducible", {

  set.seed(123)

  ## -----------------------------
  ## Minimal synthetic observed data
  ## -----------------------------

  n_years_obs <- 5
  days_per_year <- 365
  n_days_obs <- n_years_obs * days_per_year

  dates_obs <- seq.Date(
    from = as.Date("2000-01-01"),
    by = "day",
    length.out = n_days_obs
  )

  dates.d <- data.frame(
    date  = dates_obs,
    year  = rep(2000:(2000 + n_years_obs - 1), each = days_per_year),
    month = as.integer(format(dates_obs, "%m")),
    day   = as.integer(format(dates_obs, "%d")),
    wyear = rep(2000:(2000 + n_years_obs - 1), each = days_per_year)
  )

  ## Observed weather
  PRCP <- rgamma(n_days_obs, shape = 2, scale = 2)
  TEMP <- rnorm(n_days_obs, mean = 15, sd = 5)

  ## Annual observed precipitation
  ANNUAL_PRCP <- tapply(PRCP, dates.d$wyear, sum)

  ## -----------------------------
  ## Simulation setup
  ## -----------------------------

  ymax <- 2
  sim_days <- ymax * days_per_year

  sim_dates <- seq.Date(
    from = as.Date("2010-01-01"),
    by = "day",
    length.out = sim_days
  )

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = rep(2010:(2010 + ymax - 1), each = days_per_year)
  )

  ## Synthetic annual totals (target)
  PRCP_FINAL_ANNUAL_SIM <- rep(mean(ANNUAL_PRCP), ymax)

  ## -----------------------------
  ## Run resampleDates twice
  ## -----------------------------

  out1 <- resampleDates(
    PRCP_FINAL_ANNUAL_SIM = PRCP_FINAL_ANNUAL_SIM,
    ANNUAL_PRCP = ANNUAL_PRCP,
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2010,
    k1 = 1,
    ymax = ymax,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    month.start = 1,
    knn.annual.sample.num = 3,
    wet.quantile = 0.2,
    extreme.quantile = 0.8,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    alpha = 1,
    seed = 42
  )

  out2 <- resampleDates(
    PRCP_FINAL_ANNUAL_SIM = PRCP_FINAL_ANNUAL_SIM,
    ANNUAL_PRCP = ANNUAL_PRCP,
    PRCP = PRCP,
    TEMP = TEMP,
    START_YEAR_SIM = 2010,
    k1 = 1,
    ymax = ymax,
    dates.d = dates.d,
    sim.dates.d = sim.dates.d,
    month.start = 1,
    knn.annual.sample.num = 3,
    wet.quantile = 0.2,
    extreme.quantile = 0.8,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    alpha = 1,
    seed = 42
  )

  ## -----------------------------
  ## Assertions
  ## -----------------------------

  # Correct length
  expect_length(out1, sim_days)

  # Class
  expect_s3_class(out1, "Date")

  # No missing values
  expect_false(anyNA(out1))

  # All dates must come from observed record
  expect_true(all(out1 %in% dates.d$date))

  # Reproducibility
  expect_identical(out1, out2)

  # Non-degenerate behavior:
  # should not return the same date repeated
  expect_gt(length(unique(out1)), 50)
})

test_that("resampleDates performance does not regress", {

  skip_on_cran()
  skip_if(Sys.getenv("CI") == "true")  # optional but recommended

  set.seed(123)

  ## -----------------------------
  ## Moderate but realistic setup
  ## -----------------------------

  n_years_obs <- 10
  days_per_year <- 365
  n_days_obs <- n_years_obs * days_per_year

  dates_obs <- seq.Date(
    from = as.Date("1990-01-01"),
    by = "day",
    length.out = n_days_obs
  )

  dates.d <- data.frame(
    date  = dates_obs,
    year  = rep(1990:(1990 + n_years_obs - 1), each = days_per_year),
    month = as.integer(format(dates_obs, "%m")),
    day   = as.integer(format(dates_obs, "%d")),
    wyear = rep(1990:(1990 + n_years_obs - 1), each = days_per_year)
  )

  PRCP <- rgamma(n_days_obs, shape = 2, scale = 2)
  TEMP <- rnorm(n_days_obs, mean = 15, sd = 5)

  ANNUAL_PRCP <- tapply(PRCP, dates.d$wyear, sum)

  ymax <- 3
  sim_days <- ymax * days_per_year

  sim_dates <- seq.Date(
    from = as.Date("2010-01-01"),
    by = "day",
    length.out = sim_days
  )

  sim.dates.d <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = rep(2010:(2010 + ymax - 1), each = days_per_year)
  )

  PRCP_FINAL_ANNUAL_SIM <- rep(mean(ANNUAL_PRCP), ymax)

  ## -----------------------------
  ## Warm-up (important)
  ## -----------------------------

  resampleDates(
    PRCP_FINAL_ANNUAL_SIM,
    ANNUAL_PRCP,
    PRCP,
    TEMP,
    START_YEAR_SIM = 2010,
    k1 = 1,
    ymax = ymax,
    dates.d,
    sim.dates.d,
    seed = 1
  )

  ## -----------------------------
  ## Timed run
  ## -----------------------------

  t0 <- proc.time()[["elapsed"]]

  resampleDates(
    PRCP_FINAL_ANNUAL_SIM,
    ANNUAL_PRCP,
    PRCP,
    TEMP,
    START_YEAR_SIM = 2010,
    k1 = 1,
    ymax = ymax,
    dates.d,
    sim.dates.d,
    seed = 1
  )

  elapsed <- proc.time()[["elapsed"]] - t0

  ## -----------------------------
  ## Performance assertion
  ## -----------------------------

  # Generous upper bound (adjust if needed)
  # Should run well below this on normal machines
  expect_lt(elapsed, 2.0)
})
