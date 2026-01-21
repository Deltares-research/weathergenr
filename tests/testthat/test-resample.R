# Functions tested (relative paths):
# - R/resample.R: resample_weather_dates(), knn_sample(), expand_indices(),
#   estimate_monthly_markov_probs(), normalize_probs(), markov_next_state(),
#   match_transition_positions(), get_result_index()

# ---- resample_weather_dates -------------------------------------------------

testthat::test_that("resample_weather_dates returns Date vector of correct length", {
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
  obs_daily_temp <- rnorm(length(obs_dates), mean = 10, sd = 3)

  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    sim_annual_precip = 500,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    annual_knn_n      = 3,
    year_start_month  = 1,
    seed              = 123
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
  obs_daily_temp <- seq_along(obs_dates)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    sim_annual_precip = 100,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    year_start_month  = 1,
    seed              = 999
  )

  obs_year <- obs_dates_df$wyear[match(out, obs_dates_df$date)]
  obs_mon  <- as.integer(format(out, "%m"))
  obs_day  <- as.integer(format(out, "%d"))


  testthat::expect_s3_class(out, "Date")

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
  obs_daily_temp <- rnorm(length(obs_dates), 10, 3)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-10-01"), as.Date("2021-09-30"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2021
  )

  out <- resample_weather_dates(
    sim_annual_precip = 600,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
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

  obs_dates_df <- data.frame(
    date  = dates_obs,
    month = as.integer(format(dates_obs, "%m")),
    day   = as.integer(format(dates_obs, "%d")),
    wyear = as.integer(format(dates_obs, "%Y"))
  )

  obs_daily_precip <- rgamma(n, 2, 2)
  obs_daily_temp <- rnorm(n, 15, 5)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2010-01-01"), by = "day", length.out = 365)
  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2010
  )

  out1 <- resample_weather_dates(
    sim_annual_precip = mean(obs_annual_precip),
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2010,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    seed              = 42
  )

  out2 <- resample_weather_dates(
    sim_annual_precip = mean(obs_annual_precip),
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2010,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    seed              = 42
  )

  testthat::expect_identical(out1, out2)
})

testthat::test_that("resample_weather_dates rejects invalid year_start_month", {
  testthat::expect_error(
    resample_weather_dates(
      sim_annual_precip = 100,
      obs_annual_precip = 100,
      obs_daily_precip  = 1,
      obs_daily_temp    = 1,
      year_start        = 2000,
      n_years           = 1,
      obs_dates_df      = data.frame(date = Sys.Date(), month = 1, day = 1, wyear = 2000),
      sim_dates_df      = data.frame(month = 1, day = 1, wyear = 2000),
      year_start_month  = 13
    ),
    "year_start_month must be a single integer between 1 and 12"
  )
})

testthat::test_that("multi-year simulation maintains year boundaries", {
  set.seed(100)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2009-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- rgamma(length(obs_dates), shape = 2, scale = 2)
  obs_daily_temp <- rnorm(length(obs_dates), mean = 10, sd = 3)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  n_years <- 3
  sim_dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = n_years * 365)

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = rep(2020:(2020 + n_years - 1), each = 365)
  )

  sim_annual_precip <- rep(mean(obs_annual_precip), n_years)

  out <- resample_weather_dates(
    sim_annual_precip = sim_annual_precip,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = n_years,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    year_start_month  = 1,
    seed              = 456
  )

  testthat::expect_length(out, n_years * 365)
  testthat::expect_false(anyNA(out))

  sim_years <- rep(2020:(2020 + n_years - 1), each = 365)
  testthat::expect_equal(length(unique(sim_years)), n_years)
})


testthat::test_that("different year_start_month values work correctly", {
  for (start_month in c(4, 7)) {
    obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
    obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

    obs_dates_df <- data.frame(
      date  = obs_dates,
      month = as.integer(format(obs_dates, "%m")),
      day   = as.integer(format(obs_dates, "%d")),
      wyear = ifelse(
        as.integer(format(obs_dates, "%m")) >= start_month,
        as.integer(format(obs_dates, "%Y")) + 1,
        as.integer(format(obs_dates, "%Y"))
      )
    )

    obs_daily_precip <- rgamma(length(obs_dates), 2, 2)
    obs_daily_temp <- rnorm(length(obs_dates), 10, 3)
    obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

    sim_start_date <- as.Date(sprintf("2020-%02d-01", start_month))
    sim_dates <- seq.Date(sim_start_date, by = "day", length.out = 365)

    sim_dates_df <- data.frame(
      month = as.integer(format(sim_dates, "%m")),
      day   = as.integer(format(sim_dates, "%d")),
      wyear = 2021
    )

    out <- resample_weather_dates(
      sim_annual_precip = 500,
      obs_annual_precip = obs_annual_precip,
      obs_daily_precip  = obs_daily_precip,
      obs_daily_temp    = obs_daily_temp,
      year_start        = 2020,
      n_years           = 1,
      obs_dates_df      = obs_dates_df,
      sim_dates_df      = sim_dates_df,
      year_start_month  = start_month,
      seed              = 789
    )

    testthat::expect_length(out, 365)
    testthat::expect_false(anyNA(out))
    testthat::expect_equal(as.integer(format(sim_dates[1], "%m")), start_month)
  }
})

testthat::test_that("dry_spell_factor and wet_spell_factor affect results", {
  set.seed(300)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- rgamma(length(obs_dates), 2, 2)
  obs_daily_temp <- rnorm(length(obs_dates), 10, 3)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out_default <- resample_weather_dates(
    sim_annual_precip = 500,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    dry_spell_factor  = rep(1, 12),
    wet_spell_factor  = rep(1, 12),
    seed              = 111
  )

  out_dry <- resample_weather_dates(
    sim_annual_precip = 500,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    dry_spell_factor  = rep(2, 12),
    wet_spell_factor  = rep(1, 12),
    seed              = 111
  )

  testthat::expect_false(identical(out_default, out_dry))
})

testthat::test_that("wet_q and extreme_q parameters work correctly", {
  set.seed(400)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2005-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- rgamma(length(obs_dates), 2, 2)
  obs_daily_temp <- rnorm(length(obs_dates), 10, 3)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out1 <- resample_weather_dates(
    sim_annual_precip = 500,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    wet_q             = 0.1,
    extreme_q         = 0.9,
    seed              = 222
  )

  out2 <- resample_weather_dates(
    sim_annual_precip = 500,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    wet_q             = 0.3,
    extreme_q         = 0.7,
    seed              = 222
  )

  testthat::expect_length(out1, 365)
  testthat::expect_length(out2, 365)
  testthat::expect_false(anyNA(out1))
  testthat::expect_false(anyNA(out2))
})

testthat::test_that("annual_knn_n parameter works correctly", {
  set.seed(500)

  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2010-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- rgamma(length(obs_dates), 2, 2)
  obs_daily_temp <- rnorm(length(obs_dates), 10, 3)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  for (n_samples in c(3, 10, 50)) {
    out <- resample_weather_dates(
      sim_annual_precip = 500,
      obs_annual_precip = obs_annual_precip,
      obs_daily_precip  = obs_daily_precip,
      obs_daily_temp    = obs_daily_temp,
      year_start        = 2020,
      n_years           = 1,
      obs_dates_df      = obs_dates_df,
      sim_dates_df      = sim_dates_df,
      annual_knn_n      = n_samples,
      seed              = 333
    )

    testthat::expect_length(out, 365)
    testthat::expect_false(anyNA(out))
  }
})

testthat::test_that("function handles minimal observed data", {
  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2001-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- rgamma(length(obs_dates), 2, 2)
  obs_daily_temp <- rnorm(length(obs_dates), 10, 3)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  sim_dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  sim_dates <- sim_dates[format(sim_dates, "%m-%d") != "02-29"]

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = 2020
  )

  out <- resample_weather_dates(
    sim_annual_precip = 500,
    obs_annual_precip = obs_annual_precip,
    obs_daily_precip  = obs_daily_precip,
    obs_daily_temp    = obs_daily_temp,
    year_start        = 2020,
    n_years           = 1,
    obs_dates_df      = obs_dates_df,
    sim_dates_df      = sim_dates_df,
    annual_knn_n      = 2,
    seed              = 444
  )

  testthat::expect_length(out, 365)
  testthat::expect_false(anyNA(out))
})

testthat::test_that("function handles long simulations efficiently", {
  testthat::skip_on_cran()

  obs_dates <- seq.Date(as.Date("1980-01-01"), as.Date("2009-12-31"), by = "day")
  obs_dates <- obs_dates[format(obs_dates, "%m-%d") != "02-29"]

  obs_dates_df <- data.frame(
    date  = obs_dates,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d")),
    wyear = as.integer(format(obs_dates, "%Y"))
  )

  obs_daily_precip <- rgamma(length(obs_dates), 2, 2)
  obs_daily_temp <- rnorm(length(obs_dates), 10, 3)
  obs_annual_precip <- tapply(obs_daily_precip, obs_dates_df$wyear, sum)

  n_years <- 10
  sim_dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = n_years * 365)

  sim_dates_df <- data.frame(
    month = as.integer(format(sim_dates, "%m")),
    day   = as.integer(format(sim_dates, "%d")),
    wyear = rep(2020:(2020 + n_years - 1), each = 365)
  )

  sim_annual_precip <- rep(mean(obs_annual_precip), n_years)

  timing <- system.time({
    out <- resample_weather_dates(
      sim_annual_precip = sim_annual_precip,
      obs_annual_precip = obs_annual_precip,
      obs_daily_precip  = obs_daily_precip,
      obs_daily_temp    = obs_daily_temp,
      year_start        = 2020,
      n_years           = n_years,
      obs_dates_df      = obs_dates_df,
      sim_dates_df      = sim_dates_df,
      seed              = 555
    )
  })

  testthat::expect_length(out, n_years * 365)
  testthat::expect_false(anyNA(out))
  testthat::expect_lt(unname(timing["elapsed"]), 10)
})

# ---- knn_sample -------------------------------------------------------------

testthat::test_that("knn_sample returns valid indices and is reproducible with seed", {
  candidates <- matrix(1:20, ncol = 2)
  target <- c(5, 5)

  idx1 <- knn_sample(candidates, target, k = 5, n = 3, seed = 10)
  idx2 <- knn_sample(candidates, target, k = 5, n = 3, seed = 10)

  testthat::expect_length(idx1, 3)
  testthat::expect_true(all(idx1 >= 1 & idx1 <= nrow(candidates)))
  testthat::expect_identical(idx1, idx2)
})

testthat::test_that("knn_sample supports probability weighting modes", {
  set.seed(1)
  candidates <- matrix(rnorm(30), ncol = 3)
  target <- c(0, 0, 0)

  idx_rank <- knn_sample(candidates, target, k = 6, n = 4, prob = TRUE, sampling = "rank")
  idx_dist <- knn_sample(candidates, target, k = 6, n = 4, prob = TRUE, sampling = "distance")

  testthat::expect_length(idx_rank, 4)
  testthat::expect_length(idx_dist, 4)
  testthat::expect_true(all(idx_rank >= 1 & idx_rank <= nrow(candidates)))
  testthat::expect_true(all(idx_dist >= 1 & idx_dist <= nrow(candidates)))
})

# ---- expand_indices ---------------------------------------------------------

testthat::test_that("expand_indices applies offsets and respects bounds", {
  base_idx <- c(2L, 5L)
  offsets <- -1:1
  out <- weathergenr:::expand_indices(base_idx, offsets, n_max = 7L)

  testthat::expect_true(all(out > 0))
  testthat::expect_true(all((out + 1L) <= 7L))
  testthat::expect_true(all(c(1L, 2L, 3L, 4L, 5L, 6L) %in% out))
})

# ---- estimate_monthly_markov_probs -----------------------------------------

testthat::test_that("estimate_monthly_markov_probs returns normalized probabilities", {
  set.seed(2)
  n <- 100

  precip_lag1 <- rgamma(n, 2, 2)
  precip_lag0 <- rgamma(n, 2, 2)
  month_lag1 <- sample(1:12, n, replace = TRUE)
  month_lag0 <- month_lag1
  year_lag1 <- rep(2001, n)
  year_lag0 <- rep(2001, n)

  wet_threshold <- rep(quantile(precip_lag1, 0.2), 12)
  extreme_threshold <- rep(quantile(precip_lag1, 0.8), 12)

  sim_month <- rep(1:12, length.out = 365)
  sim_wyear <- rep(2001, 365)

  out <- estimate_monthly_markov_probs(
    precip_lag0 = precip_lag0,
    precip_lag1 = precip_lag1,
    month_lag0 = month_lag0,
    month_lag1 = month_lag1,
    year_lag0 = year_lag0,
    year_lag1 = year_lag1,
    wet_threshold = wet_threshold,
    extreme_threshold = extreme_threshold,
    month_order = 1:12,
    sim_month = sim_month,
    sim_wyear = sim_wyear,
    year_idx = 1,
    sim_start_year = 2001,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = 365,
    dirichlet_alpha = 1
  )

  testthat::expect_true(all(c("p00_final", "p01_final", "p02_final") %in% names(out)))
  testthat::expect_length(out$p00_final, 365)
  testthat::expect_true(all(out$p00_final >= 0 & out$p00_final <= 1))

  row_sum <- out$p00_final + out$p01_final + out$p02_final
  testthat::expect_true(all(abs(row_sum - 1) < 1e-8 | is.na(row_sum)))
})

# ---- normalize_probs --------------------------------------------------------

testthat::test_that("normalize_probs handles non-finite values and zero mass", {
  out <- weathergenr:::normalize_probs(c(0.2, NA, -1, 0.3))
  testthat::expect_equal(sum(out), 1)
  testthat::expect_true(all(out >= 0))

  out_zero <- weathergenr:::normalize_probs(c(0, 0, 0))
  testthat::expect_equal(sum(out_zero), 1)
  testthat::expect_true(all(out_zero == rep(1 / 3, 3)))
})

# ---- markov_next_state ------------------------------------------------------

testthat::test_that("markov_next_state returns valid state and respects bounds", {
  p00 <- rep(0.7, 5)
  p01 <- rep(0.2, 5)
  p10 <- rep(0.3, 5)
  p11 <- rep(0.4, 5)
  p20 <- rep(0.1, 5)
  p21 <- rep(0.3, 5)

  out_low <- markov_next_state(0L, u_rand = -0.2, idx = 1, p00, p01, p10, p11, p20, p21)
  out_high <- markov_next_state(1L, u_rand = 1.2, idx = 10, p00, p01, p10, p11, p20, p21)

  testthat::expect_true(out_low %in% 0:2)
  testthat::expect_true(out_high %in% 0:2)
})

# ---- match_transition_positions ---------------------------------------------

testthat::test_that("match_transition_positions identifies expected transitions", {
  precip_vec <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
  day0_idx <- 1:(length(precip_vec) - 1)
  wet_threshold <- 1
  extreme_threshold <- 20

  idx_0_1 <- match_transition_positions(0, 1, precip_vec, day0_idx, wet_threshold, extreme_threshold)
  idx_1_2 <- match_transition_positions(1, 2, precip_vec, day0_idx, wet_threshold, extreme_threshold)

  testthat::expect_true(length(idx_0_1) > 0)
  testthat::expect_true(length(idx_1_2) > 0)
})

# ---- get_result_index -------------------------------------------------------

testthat::test_that("get_result_index returns provided index or samples", {
  candidate_precip <- c(0, 5, 10, 20)

  testthat::expect_identical(get_result_index(2, candidate_precip), 2L)

  set.seed(1)
  idx <- get_result_index(10, candidate_precip)
  testthat::expect_true(idx %in% seq_along(candidate_precip))

  testthat::expect_true(is.na(get_result_index(1, numeric(0))))
})
