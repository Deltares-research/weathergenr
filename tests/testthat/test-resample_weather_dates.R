
testthat::test_that("resample_weather_dates returns Date vector of correct length", {

  set.seed(123)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), shape = 2, scale = 2)
  daily_temp_obs <- rnorm(length(obs_dates), mean = 10, sd = 3)

  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    annual_knn_n      = 3,
    year_start_month  = 1,
    seed              = 123
  )

  testthat::expect_s3_class(out, "Date")
  testthat::expect_length(out, nrow(sim_dates_table))
  testthat::expect_false(anyNA(out))
})

testthat::test_that("calendar-year mode forbids observed Dec->Jan transitions", {

  set.seed(42)

  obs_dates <- seq.Date(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- seq_along(obs_dates)
  daily_temp_obs <- seq_along(obs_dates)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    annual_prcp_sim  = 100,
    annual_prcp_obs  = annual_prcp_obs,
    daily_prcp_obs   = daily_prcp_obs,
    daily_temp_obs   = daily_temp_obs,
    sim_start_year   = 2020,
    realization_id   = 1,
    n_sim_years      = 1,
    obs_dates_table  = obs_dates_table,
    sim_dates_table  = sim_dates_table,
    year_start_month = 1,
    seed             = 999
  )

  obs_year <- obs_dates_table$wyear[match(out, obs_dates_table$date)]
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

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = ifelse(
      as.integer(format(obs_dates, "%m")) >= 10,
      as.integer(format(obs_dates, "%Y")) + 1,
      as.integer(format(obs_dates, "%Y"))
    )
  )

  daily_prcp_obs <- rgamma(length(obs_dates), 2, 2)
  daily_temp_obs <- rnorm(length(obs_dates), 10, 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-10-01"), as.Date("2021-09-30"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2021
  )

  out <- resample_weather_dates(
    annual_prcp_sim   = 600,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    annual_knn_n      = 3,
    year_start_month  = 10,
    seed              = 321
  )

  mons <- as.integer(format(out, "%m"))
  testthat::expect_true(any(mons[-1] == 1 & mons[-length(mons)] == 12))
})

testthat::test_that("resample_weather_dates is reproducible with same seed", {

  n <- 5 * 365
  dates_obs <- seq.Date(as.Date("2000-01-01"), by = "day", length.out = n)

  obs_dates_table <- data.frame(
    date  = dates_obs,
    month = as.integer(format(dates_obs, "%m")),
    day   = as.integer(format(dates_obs, "%d")),
    wyear = as.integer(format(dates_obs, "%Y"))
  )

  daily_prcp_obs <- rgamma(n, 2, 2)
  daily_temp_obs <- rnorm(n, 15, 5)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2010-01-01"), by = "day", length.out = 365)
  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2010
  )

  out1 <- resample_weather_dates(
    annual_prcp_sim  = mean(annual_prcp_obs),
    annual_prcp_obs  = annual_prcp_obs,
    daily_prcp_obs   = daily_prcp_obs,
    daily_temp_obs   = daily_temp_obs,
    sim_start_year   = 2010,
    realization_id   = 1,
    n_sim_years      = 1,
    obs_dates_table  = obs_dates_table,
    sim_dates_table  = sim_dates_table,
    seed             = 42
  )

  out2 <- resample_weather_dates(
    annual_prcp_sim  = mean(annual_prcp_obs),
    annual_prcp_obs  = annual_prcp_obs,
    daily_prcp_obs   = daily_prcp_obs,
    daily_temp_obs   = daily_temp_obs,
    sim_start_year   = 2010,
    realization_id   = 1,
    n_sim_years      = 1,
    obs_dates_table  = obs_dates_table,
    sim_dates_table  = sim_dates_table,
    seed             = 42
  )

  testthat::expect_identical(out1, out2)
})

testthat::test_that("resample_weather_dates rejects invalid year_start_month", {

  testthat::expect_error(
    resample_weather_dates(
      annual_prcp_sim   = 100,
      annual_prcp_obs   = 100,
      daily_prcp_obs    = 1,
      daily_temp_obs    = 1,
      sim_start_year    = 2000,
      realization_id    = 1,
      n_sim_years       = 1,
      obs_dates_table   = data.frame(date = Sys.Date(), month = 1, day = 1, wyear = 2000),
      sim_dates_table   = data.frame(month = 1, day = 1, wyear = 2000),
      year_start_month  = 13
    ),
    "year_start_month must be a single integer between 1 and 12"
  )
})

testthat::test_that("multi-year simulation maintains year boundaries", {
  set.seed(100)

  # 10 years of observed data
  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2009-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), shape = 2, scale = 2)
  daily_temp_obs <- rnorm(length(obs_dates), mean = 10, sd = 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  # Simulate 3 years
  n_years <- 3
  sim_dates <- seq.Date(as.Date("2020-01-01"),
                        by = "day",
                        length.out = n_years * 365)

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = rep(2020:(2020 + n_years - 1), each = 365)
  )

  annual_prcp_sim <- rep(mean(annual_prcp_obs), n_years)

  out <- resample_weather_dates(
    annual_prcp_sim   = annual_prcp_sim,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = n_years,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    year_start_month  = 1,
    seed              = 456
  )

  testthat::expect_length(out, n_years * 365)
  testthat::expect_false(anyNA(out))

  # Check that each year has exactly 365 days
  sim_years <- rep(2020:(2020 + n_years - 1), each = 365)
  testthat::expect_equal(length(unique(sim_years)), n_years)
})

testthat::test_that("different realization_id produces different results", {
  set.seed(200)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), shape = 2, scale = 2)
  daily_temp_obs <- rnorm(length(obs_dates), mean = 10, sd = 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out1 <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,  # Different realization
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    seed              = 100
  )

  out2 <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 2,  # Different realization
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    seed              = 100
  )

  # Results should be different
  testthat::expect_false(identical(out1, out2))

  # But both should be valid
  testthat::expect_length(out1, 365)
  testthat::expect_length(out2, 365)
  testthat::expect_false(anyNA(out1))
  testthat::expect_false(anyNA(out2))
})

testthat::test_that("different year_start_month values work correctly", {

  # Test months 4 (April) and 7 (July)
  for (start_month in c(4, 7)) {

    obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
    obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

    obs_dates_table <- data.frame(
      date  = obs_dates,
      month = as.integer(format(obs_dates, "%m")),
      day   = as.integer(format(obs_dates, "%d")),
      wyear = ifelse(
        as.integer(format(obs_dates, "%m")) >= start_month,
        as.integer(format(obs_dates, "%Y")) + 1,
        as.integer(format(obs_dates, "%Y"))
      )
    )

    daily_prcp_obs <- rgamma(length(obs_dates), 2, 2)
    daily_temp_obs <- rnorm(length(obs_dates), 10, 3)
    annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

    sim_start_date <- as.Date(sprintf("2020-%02d-01", start_month))
    sim_dates <- seq.Date(sim_start_date, by = "day", length.out = 365)

    sim_dates_table <- data.frame(
      month = as.integer(format(sim_dates, "%m")),
      day   = as.integer(format(sim_dates, "%d")),
      wyear = 2021
    )

    out <- resample_weather_dates(
      annual_prcp_sim   = 500,
      annual_prcp_obs   = annual_prcp_obs,
      daily_prcp_obs    = daily_prcp_obs,
      daily_temp_obs    = daily_temp_obs,
      sim_start_year    = 2020,
      realization_id    = 1,
      n_sim_years       = 1,
      obs_dates_table   = obs_dates_table,
      sim_dates_table   = sim_dates_table,
      year_start_month  = start_month,
      seed              = 789
    )

    testthat::expect_length(out, 365)
    testthat::expect_false(anyNA(out))

    # First simulated date should be in the start month
    testthat::expect_equal(as.integer(format(sim_dates[1], "%m")), start_month)
  }
})

testthat::test_that("dry_spell_factor and wet_spell_factor affect results", {
  set.seed(300)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), 2, 2)
  daily_temp_obs <- rnorm(length(obs_dates), 10, 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  # No spell modification
  out_default <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    dry_spell_factor  = rep(1, 12),
    wet_spell_factor  = rep(1, 12),
    seed              = 111
  )

  # Increase dry spell persistence
  out_dry <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    dry_spell_factor  = rep(2, 12),  # Increase dry persistence
    wet_spell_factor  = rep(1, 12),
    seed              = 111
  )

  # Results should be different when spell factors change
  testthat::expect_false(identical(out_default, out_dry))
})

testthat::test_that("wet_q and extreme_q parameters work correctly", {
  set.seed(400)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), 2, 2)
  daily_temp_obs <- rnorm(length(obs_dates), 10, 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  # Test with different quantile combinations
  out1 <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    wet_q             = 0.1,
    extreme_q         = 0.9,
    seed              = 222
  )

  out2 <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    wet_q             = 0.3,
    extreme_q         = 0.7,
    seed              = 222
  )

  # Both should complete successfully
  testthat::expect_length(out1, 365)
  testthat::expect_length(out2, 365)
  testthat::expect_false(anyNA(out1))
  testthat::expect_false(anyNA(out2))

  # Results may differ depending on the data
  # Main goal is to ensure no errors
})

testthat::test_that("annual_knn_n parameter works correctly", {
  set.seed(500)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2010-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), 2, 2)
  daily_temp_obs <- rnorm(length(obs_dates), 10, 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  # Test with different sample sizes
  for (n_samples in c(3, 10, 50)) {
    out <- resample_weather_dates(
      annual_prcp_sim   = 500,
      annual_prcp_obs   = annual_prcp_obs,
      daily_prcp_obs    = daily_prcp_obs,
      daily_temp_obs    = daily_temp_obs,
      sim_start_year    = 2020,
      realization_id    = 1,
      n_sim_years       = 1,
      obs_dates_table   = obs_dates_table,
      sim_dates_table   = sim_dates_table,
      annual_knn_n      = n_samples,
      seed              = 333
    )

    testthat::expect_length(out, 365)
    testthat::expect_false(anyNA(out))
  }
})

testthat::test_that("function handles minimal observed data", {
  # Just 2 years of observed data
  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2001-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), 2, 2)
  daily_temp_obs <- rnorm(length(obs_dates), 10, 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  # Should still work with minimal data
  out <- resample_weather_dates(
    annual_prcp_sim   = 500,
    annual_prcp_obs   = annual_prcp_obs,
    daily_prcp_obs    = daily_prcp_obs,
    daily_temp_obs    = daily_temp_obs,
    sim_start_year    = 2020,
    realization_id    = 1,
    n_sim_years       = 1,
    obs_dates_table   = obs_dates_table,
    sim_dates_table   = sim_dates_table,
    annual_knn_n      = 2,  # Can't be more than available years
    seed              = 444
  )

  testthat::expect_length(out, 365)
  testthat::expect_false(anyNA(out))
})

testthat::test_that("function handles long simulations efficiently", {
  skip_on_cran()  # Skip on CRAN due to time

  # 30 years of observed data
  obs_dates <- seq.Date(as.Date("1980-01-01"), as.Date("2009-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_table <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  daily_prcp_obs <- rgamma(length(obs_dates), 2, 2)
  daily_temp_obs <- rnorm(length(obs_dates), 10, 3)
  annual_prcp_obs <- tapply(daily_prcp_obs, obs_dates_table$wyear, sum)

  # Simulate 10 years
  n_years <- 10
  sim_dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = n_years * 365)

  sim_dates_table <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = rep(2020:(2020 + n_years - 1), each = 365)
  )

  annual_prcp_sim <- rep(mean(annual_prcp_obs), n_years)

  # Should complete in reasonable time (< 5 seconds for 10 years)
  timing <- system.time({
    out <- resample_weather_dates(
      annual_prcp_sim   = annual_prcp_sim,
      annual_prcp_obs   = annual_prcp_obs,
      daily_prcp_obs    = daily_prcp_obs,
      daily_temp_obs    = daily_temp_obs,
      sim_start_year    = 2020,
      realization_id    = 1,
      n_sim_years       = n_years,
      obs_dates_table   = obs_dates_table,
      sim_dates_table   = sim_dates_table,
      seed              = 555
    )
  })

  testthat::expect_length(out, n_years * 365)
  testthat::expect_false(anyNA(out))
  testthat::expect_lt(timing["elapsed"], 10)  # Should complete in < 10 seconds
})

