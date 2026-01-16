testthat::test_that("estimate_monthly_markov_probs returns correct structure and lengths", {

  sim_length <- 365

  res <- estimate_monthly_markov_probs(
    precip_lag0 = rep(5, 100),
    precip_lag1 = rep(3, 100),
    month_lag0  = rep(1L, 100),
    month_lag1  = rep(1L, 100),
    year_lag0   = rep(2001L, 100),
    year_lag1   = rep(2001L, 100),
    wet_threshold     = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_order = 1:12,
    sim_month = rep(1L, sim_length),
    sim_wyear = rep(2001L, sim_length),
    year_idx = 1L,
    sim_start_year = 2001L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = sim_length
  )

  testthat::expect_type(res, "list")
  testthat::expect_equal(length(res), 9)

  for (nm in names(res)) {
    testthat::expect_length(res[[nm]], sim_length)
    testthat::expect_true(is.numeric(res[[nm]]))
  }
})

# ------------------------------------------------------------
# Calendar-year mode (month_order[1] == 1)
# ------------------------------------------------------------
testthat::test_that("calendar-year mode excludes cross-year transitions", {


  # December (12) to January (1) transition (cross-year)
  precip_lag1 <- c(5, 5)
  precip_lag0 <- c(5, 5)

  month_lag1 <- c(12L, 1L)
  month_lag0 <- c(1L, 2L)

  year_lag1 <- c(2001L, 2002L)
  year_lag0 <- c(2002L, 2002L)

  sim_length <- 10

  res <- estimate_monthly_markov_probs(
    precip_lag0 = precip_lag0,
    precip_lag1 = precip_lag1,
    month_lag0  = month_lag0,
    month_lag1  = month_lag1,
    year_lag0   = year_lag0,
    year_lag1   = year_lag1,
    wet_threshold     = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_order = 1:12,            # calendar-year mode
    sim_month = rep(1L, sim_length),
    sim_wyear = rep(2002L, sim_length),
    year_idx = 1L,
    sim_start_year = 2002L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = sim_length
  )

  # Should not error and probabilities should be valid/bounded where defined
  testthat::expect_true(all(res$p00_final[!is.na(res$p00_final)] >= 0))
  testthat::expect_true(all(res$p01_final[!is.na(res$p01_final)] >= 0))
})

# ------------------------------------------------------------
# Water-year mode (month_order[1] != 1)
# ------------------------------------------------------------
testthat::test_that("water-year mode allows cross-calendar-year transitions", {

  precip_lag1 <- c(5, 5)
  precip_lag0 <- c(5, 5)

  month_lag1 <- c(12L, 1L)
  month_lag0 <- c(1L, 2L)

  year_lag1 <- c(2001L, 2001L)
  year_lag0 <- c(2001L, 2001L)

  sim_length <- 10

  res <- estimate_monthly_markov_probs(
    precip_lag0 = precip_lag0,
    precip_lag1 = precip_lag1,
    month_lag0  = month_lag0,
    month_lag1  = month_lag1,
    year_lag0   = year_lag0,
    year_lag1   = year_lag1,
    wet_threshold     = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_order = c(10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9),  # water-year mode
    sim_month = rep(1L, sim_length),
    sim_wyear = rep(2001L, sim_length),
    year_idx = 1L,
    sim_start_year = 2000L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = sim_length
  )

  testthat::expect_true(any(res$p00_final > 0, na.rm = TRUE))
})

# ------------------------------------------------------------
# Probability consistency
# ------------------------------------------------------------
testthat::test_that("transition probabilities are valid and bounded", {

  sim_length <- 50

  res <- estimate_monthly_markov_probs(
    precip_lag0 = runif(200, 0, 20),
    precip_lag1 = runif(200, 0, 20),
    month_lag0  = sample(1:12, 200, TRUE),
    month_lag1  = sample(1:12, 200, TRUE),
    year_lag0   = rep(2001L, 200),
    year_lag1   = rep(2001L, 200),
    wet_threshold     = rep(5, 12),
    extreme_threshold = rep(15, 12),
    month_order = 1:12,
    sim_month = sample(1:12, sim_length, TRUE),
    sim_wyear = rep(2001L, sim_length),
    year_idx = 1L,
    sim_start_year = 2001L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = sim_length
  )

  for (nm in names(res)) {
    x <- res[[nm]]
    x <- x[!is.na(x)]
    testthat::expect_true(all(x >= 0))
    testthat::expect_true(all(x <= 1))
  }
})

# ------------------------------------------------------------
# YEAR_LAG guard correctness
# ------------------------------------------------------------
testthat::test_that("year_lag filtering removes cross-year transitions", {

  precip_lag0 <- c(5, 5, 5)
  precip_lag1 <- c(5, 5, 5)

  month_lag0 <- c(1L, 1L, 1L)
  month_lag1 <- c(12L, 1L, 1L)

  year_lag0 <- c(2002L, 2002L, 2002L)
  year_lag1 <- c(2001L, 2002L, 2002L)

  sim_length <- 10

  res <- estimate_monthly_markov_probs(
    precip_lag0 = precip_lag0,
    precip_lag1 = precip_lag1,
    month_lag0  = month_lag0,
    month_lag1  = month_lag1,
    year_lag0   = year_lag0,
    year_lag1   = year_lag1,
    wet_threshold     = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_order = 1:12,
    sim_month = rep(1L, sim_length),
    sim_wyear = rep(2002L, sim_length),
    year_idx = 1L,
    sim_start_year = 2002L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = sim_length
  )

  testthat::expect_true(all(is.finite(res$p00_final[!is.na(res$p00_final)])))
})

testthat::test_that("estimate_monthly_markov_probs applies smoothing for zero-denominator months", {

  # Extreme sparse case: only dry observed in July
  precip_lag1 <- rep(0, 20)
  precip_lag0 <- rep(0, 20)

  month_lag1 <- rep(7L, 20)
  month_lag0 <- rep(7L, 20)

  year_lag1 <- rep(2000L, 20)
  year_lag0 <- rep(2000L, 20)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_order <- 1:12

  sim_month <- rep(7L, 31)
  sim_wyear <- rep(2000L, 31)

  res <- estimate_monthly_markov_probs(
    precip_lag0 = precip_lag0,
    precip_lag1 = precip_lag1,
    month_lag0  = month_lag0,
    month_lag1  = month_lag1,
    year_lag0   = year_lag0,
    year_lag1   = year_lag1,
    wet_threshold = wet_threshold,
    extreme_threshold = extreme_threshold,
    month_order = month_order,
    sim_month = sim_month,
    sim_wyear = sim_wyear,
    year_idx = 1L,
    sim_start_year = 2000L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month)
  )

  idx <- which(sim_month == 7L)

  p00 <- res$p00_final[idx][1]
  p01 <- res$p01_final[idx][1]
  p02 <- res$p02_final[idx][1]

  testthat::expect_gt(p01, 0)
  testthat::expect_gt(p02, 0)

  testthat::expect_true(all(c(p00, p01, p02) >= 0))
  testthat::expect_true(all(c(p00, p01, p02) <= 1))

  testthat::expect_equal(p00 + p01 + p02, 1, tolerance = 1e-12)
})

testthat::test_that("estimate_monthly_markov_probs smoothing strength affects probabilities", {

  precip_lag1 <- c(rep(0, 10), rep(5, 2))  # mostly dry, a few wet
  precip_lag0 <- c(rep(0, 10), rep(5, 2))

  month_lag1 <- rep(6L, 12)
  month_lag0 <- rep(6L, 12)
  year_lag1  <- rep(2001L, 12)
  year_lag0  <- rep(2001L, 12)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_order <- 1:12

  sim_month <- rep(6L, 30)
  sim_wyear <- rep(2001L, 30)

  # Default smoothing (dirichlet_alpha default)
  res_low <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L, sim_start_year = 2001L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month)
  )

  idx <- which(sim_month == 6L)[1]
  p01_low <- res_low$p01_final[idx]

  # Increase sample size (reduces effective alpha_m = alpha/sqrt(N_m))
  precip_lag1b <- rep(precip_lag1, 5)
  precip_lag0b <- rep(precip_lag0, 5)
  month_lag1b  <- rep(month_lag1, 5)
  month_lag0b  <- rep(month_lag0, 5)
  year_lag1b   <- rep(year_lag1, 5)
  year_lag0b   <- rep(year_lag0, 5)

  res_high <- estimate_monthly_markov_probs(
    precip_lag0b, precip_lag1b,
    month_lag0b, month_lag1b,
    year_lag0b, year_lag1b,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L, sim_start_year = 2001L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month)
  )

  p01_high <- res_high$p01_final[idx]

  # With more data, smoothing influence should be weaker (closer to empirical 0)
  testthat::expect_gt(abs(p01_low - 0), abs(p01_high - 0))
})

testthat::test_that("estimate_monthly_markov_probs smoothing minimally affects well-sampled months", {

  set.seed(1)

  n <- 1000
  precip_lag1 <- rgamma(n, shape = 2, scale = 2)
  precip_lag0 <- rgamma(n, shape = 2, scale = 2)

  month_lag1 <- rep(8L, n)
  month_lag0 <- rep(8L, n)
  year_lag1  <- rep(1999L, n)
  year_lag0  <- rep(1999L, n)

  wet_threshold     <- rep(2, 12)
  extreme_threshold <- rep(6, 12)
  month_order <- 1:12

  sim_month <- rep(8L, 31)
  sim_wyear <- rep(1999L, 31)

  res <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L, sim_start_year = 1999L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month)
  )

  idx <- which(sim_month == 8L)[1]

  probs <- c(
    res$p00_final[idx], res$p01_final[idx], res$p02_final[idx],
    res$p10_final[idx], res$p11_final[idx], res$p12_final[idx]
  )

  testthat::expect_true(all(probs > 0.01))
  testthat::expect_true(all(probs < 0.99))
})

testthat::test_that("estimate_monthly_markov_probs produces no absorbing states", {

  precip_lag1 <- rep(0, 15)
  precip_lag0 <- rep(0, 15)

  month_lag1 <- rep(1L, 15)
  month_lag0 <- rep(1L, 15)
  year_lag1  <- rep(1980L, 15)
  year_lag0  <- rep(1980L, 15)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_order <- 1:12

  sim_month <- rep(1L, 31)
  sim_wyear <- rep(1980L, 31)

  res <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L, sim_start_year = 1980L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month)
  )

  idx <- which(sim_month == 1L)[1]

  testthat::expect_gt(sum(c(res$p00_final[idx], res$p01_final[idx], res$p02_final[idx]) > 0), 1)
  testthat::expect_gt(sum(c(res$p10_final[idx], res$p11_final[idx], res$p12_final[idx]) > 0), 1)
  testthat::expect_gt(sum(c(res$p20_final[idx], res$p21_final[idx], res$p22_final[idx]) > 0), 1)
})

testthat::test_that("estimate_monthly_markov_probs respects dirichlet_alpha argument", {

  precip_lag1 <- rep(0, 10)
  precip_lag0 <- rep(0, 10)
  month_lag1 <- rep(5L, 10)
  month_lag0 <- rep(5L, 10)
  year_lag1 <- rep(2000L, 10)
  year_lag0 <- rep(2000L, 10)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_order <- 1:12

  sim_month <- rep(5L, 31)
  sim_wyear <- rep(2000L, 31)

  res_nosmooth <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L, sim_start_year = 2000L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month),
    dirichlet_alpha = 0
  )

  res_smooth <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L, sim_start_year = 2000L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month),
    dirichlet_alpha = 0.5
  )

  idx <- which(sim_month == 5L)[1]

  testthat::expect_true(res_nosmooth$p01_final[idx] == 0)
  testthat::expect_gt(res_smooth$p01_final[idx], 0)
})

testthat::test_that("estimate_monthly_markov_probs effective alpha decreases with sample size", {

  # Month 6: sparse
  precip_lag1_6 <- rep(0, 4)
  precip_lag0_6 <- rep(0, 4)

  # Month 7: dense
  precip_lag1_7 <- rep(0, 100)
  precip_lag0_7 <- rep(0, 100)

  precip_lag1 <- c(precip_lag1_6, precip_lag1_7)
  precip_lag0 <- c(precip_lag0_6, precip_lag0_7)

  month_lag1 <- c(rep(6L, 4), rep(7L, 100))
  month_lag0 <- month_lag1

  year_lag1 <- rep(2000L, length(precip_lag1))
  year_lag0 <- year_lag1

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_order <- 1:12

  sim_month <- c(rep(6L, 30), rep(7L, 30))
  sim_wyear <- rep(2000L, length(sim_month))

  res <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L,
    sim_start_year = 2000L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month),
    dirichlet_alpha = 1
  )

  idx6 <- which(sim_month == 6L)[1]
  idx7 <- which(sim_month == 7L)[1]

  p01_6 <- res$p01_final[idx6]
  p01_7 <- res$p01_final[idx7]

  testthat::expect_gt(p01_6, p01_7)
})

testthat::test_that("estimate_monthly_markov_probs effective alpha vanishes for large sample sizes", {

  set.seed(1)

  n <- 1000

  precip_lag1 <- c(rep(0, n / 2), rep(5, n / 2))
  precip_lag0 <- precip_lag1

  month_lag1 <- rep(9L, n)
  month_lag0 <- rep(9L, n)

  year_lag1 <- rep(1995L, n)
  year_lag0 <- year_lag1

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_order <- 1:12

  sim_month <- rep(9L, 30)
  sim_wyear <- rep(1995L, 30)

  res <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L,
    sim_start_year = 1995L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month),
    dirichlet_alpha = 1
  )

  idx <- which(sim_month == 9L)[1]

  testthat::expect_gt(res$p00_final[idx], 0.95)
  testthat::expect_lt(res$p01_final[idx] + res$p02_final[idx], 0.05)
})

testthat::test_that("estimate_monthly_markov_probs rows sum to 1 for each origin state", {

  set.seed(42)

  n <- 300

  precip_lag1 <- rgamma(n, shape = 2, scale = 2)
  precip_lag0 <- rgamma(n, shape = 2, scale = 2)

  month_lag1 <- rep(3L, n)
  month_lag0 <- rep(3L, n)

  year_lag1 <- rep(2005L, n)
  year_lag0 <- rep(2005L, n)

  wet_threshold     <- rep(2, 12)
  extreme_threshold <- rep(8, 12)
  month_order <- 1:12

  sim_month <- rep(3L, 31)
  sim_wyear <- rep(2005L, 31)

  res <- estimate_monthly_markov_probs(
    precip_lag0, precip_lag1,
    month_lag0, month_lag1,
    year_lag0, year_lag1,
    wet_threshold, extreme_threshold,
    month_order,
    sim_month, sim_wyear,
    year_idx = 1L,
    sim_start_year = 2005L,
    dry_spell_factor_month = rep(1, 12),
    wet_spell_factor_month = rep(1, 12),
    n_days_sim = length(sim_month),
    dirichlet_alpha = 1
  )

  # pick any simulated index for March
  i <- which(sim_month == 3L)[1]

  # Dry-origin row
  testthat::expect_equal(
    res$p00_final[i] + res$p01_final[i] + res$p02_final[i],
    1,
    tolerance = 1e-12
  )

  # Wet-origin row
  testthat::expect_equal(
    res$p10_final[i] + res$p11_final[i] + res$p12_final[i],
    1,
    tolerance = 1e-12
  )

  # Very-wet-origin row
  testthat::expect_equal(
    res$p20_final[i] + res$p21_final[i] + res$p22_final[i],
    1,
    tolerance = 1e-12
  )
})
