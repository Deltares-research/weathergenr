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
