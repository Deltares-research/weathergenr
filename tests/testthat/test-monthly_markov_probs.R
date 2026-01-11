testthat::test_that("monthly_markov_probs returns correct structure and lengths", {

  SIM_LENGTH <- 365

  res <- monthly_markov_probs(
    precip.lag0 = rep(5, 100),
    precip.lag1 = rep(3, 100),
    month.lag0  = rep(1L, 100),
    month.lag1  = rep(1L, 100),
    year.lag0   = rep(2001L, 100),
    year.lag1   = rep(2001L, 100),
    wet.threshold     = rep(2, 12),
    extreme.threshold = rep(10, 12),
    month.list = 1:12,
    sim.months = rep(1L, SIM_LENGTH),
    sim.water.years = rep(2001L, SIM_LENGTH),
    year.idx = 1L,
    sim.start.year = 2001L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = SIM_LENGTH
  )

  testthat::expect_type(res, "list")
  testthat::expect_equal(length(res), 9)

  for (nm in names(res)) {
    testthat::expect_length(res[[nm]], SIM_LENGTH)
    testthat::expect_true(is.numeric(res[[nm]]))
  }
})

# ------------------------------------------------------------
# Calendar-year mode (month.list[1] == 1)
# ------------------------------------------------------------
testthat::test_that("calendar-year mode excludes cross-year transitions", {

  # December (12) to January (1) transition (cross-year)
  precip.lag1 <- c(5, 5)
  precip.lag0 <- c(5, 5)

  month.lag1 <- c(12L, 1L)
  month.lag0 <- c(1L, 2L)

  year.lag1 <- c(2001L, 2002L)
  year.lag0 <- c(2002L, 2002L)

  SIM_LENGTH <- 10

  res <- monthly_markov_probs(
    precip.lag0 = precip.lag0,
    precip.lag1 = precip.lag1,
    month.lag0  = month.lag0,
    month.lag1  = month.lag1,
    year.lag0   = year.lag0,
    year.lag1   = year.lag1,
    wet.threshold     = rep(2, 12),
    extreme.threshold = rep(10, 12),
    month.list = 1:12,            # calendar-year mode
    sim.months = rep(1L, SIM_LENGTH),
    sim.water.years = rep(2002L, SIM_LENGTH),
    year.idx = 1L,
    sim.start.year = 2002L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = SIM_LENGTH
  )

  # Should not error and probabilities should be valid/bounded where defined
  testthat::expect_true(all(res$p00_final[!is.na(res$p00_final)] >= 0))
  testthat::expect_true(all(res$p01_final[!is.na(res$p01_final)] >= 0))
})

# ------------------------------------------------------------
# Water-year mode (month.list[1] != 1)
# ------------------------------------------------------------
testthat::test_that("water-year mode allows cross-calendar-year transitions", {

  precip.lag1 <- c(5, 5)
  precip.lag0 <- c(5, 5)

  month.lag1 <- c(12L, 1L)
  month.lag0 <- c(1L, 2L)

  year.lag1 <- c(2001L, 2001L)
  year.lag0 <- c(2001L, 2001L)

  SIM_LENGTH <- 10

  res <- monthly_markov_probs(
    precip.lag0 = precip.lag0,
    precip.lag1 = precip.lag1,
    month.lag0  = month.lag0,
    month.lag1  = month.lag1,
    year.lag0   = year.lag0,
    year.lag1   = year.lag1,
    wet.threshold     = rep(2, 12),
    extreme.threshold = rep(10, 12),
    month.list = c(10,11,12,1,2,3,4,5,6,7,8,9),  # water-year mode
    sim.months = rep(1L, SIM_LENGTH),
    sim.water.years = rep(2001L, SIM_LENGTH),
    year.idx = 1L,
    sim.start.year = 2000L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = SIM_LENGTH
  )

  testthat::expect_true(any(res$p00_final > 0, na.rm = TRUE))
})

# ------------------------------------------------------------
# Probability consistency
# ------------------------------------------------------------
testthat::test_that("transition probabilities are valid and bounded", {

  SIM_LENGTH <- 50

  res <- monthly_markov_probs(
    precip.lag0 = runif(200, 0, 20),
    precip.lag1 = runif(200, 0, 20),
    month.lag0  = sample(1:12, 200, TRUE),
    month.lag1  = sample(1:12, 200, TRUE),
    year.lag0   = rep(2001L, 200),
    year.lag1   = rep(2001L, 200),
    wet.threshold     = rep(5, 12),
    extreme.threshold = rep(15, 12),
    month.list = 1:12,
    sim.months = sample(1:12, SIM_LENGTH, TRUE),
    sim.water.years = rep(2001L, SIM_LENGTH),
    year.idx = 1L,
    sim.start.year = 2001L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = SIM_LENGTH
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
testthat::test_that("YEAR_LAG filtering removes cross-year transitions", {

  precip.lag0 <- c(5, 5, 5)
  precip.lag1 <- c(5, 5, 5)

  month.lag0 <- c(1L, 1L, 1L)
  month.lag1 <- c(12L, 1L, 1L)

  year.lag0 <- c(2002L, 2002L, 2002L)
  year.lag1 <- c(2001L, 2002L, 2002L)

  SIM_LENGTH <- 10

  res <- monthly_markov_probs(
    precip.lag0 = precip.lag0,
    precip.lag1 = precip.lag1,
    month.lag0  = month.lag0,
    month.lag1  = month.lag1,
    year.lag0   = year.lag0,
    year.lag1   = year.lag1,
    wet.threshold     = rep(2, 12),
    extreme.threshold = rep(10, 12),
    month.list = 1:12,
    sim.months = rep(1L, SIM_LENGTH),
    sim.water.years = rep(2002L, SIM_LENGTH),
    year.idx = 1L,
    sim.start.year = 2002L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = SIM_LENGTH
  )

  testthat::expect_true(all(is.finite(res$p00_final[!is.na(res$p00_final)])))
})

testthat::test_that("monthly_markov_probs applies smoothing for zero-denominator months", {

  # Extreme sparse case: only dry observed in July
  precip.lag1 <- rep(0, 20)
  precip.lag0 <- rep(0, 20)

  month.lag1 <- rep(7L, 20)
  month.lag0 <- rep(7L, 20)

  year.lag1 <- rep(2000L, 20)
  year.lag0 <- rep(2000L, 20)

  wet.threshold     <- rep(1, 12)
  extreme.threshold <- rep(10, 12)
  month.list <- 1:12

  sim.months <- rep(7L, 31)
  sim.water.years <- rep(2000L, 31)

  res <- monthly_markov_probs(
    precip.lag0 = precip.lag0,
    precip.lag1 = precip.lag1,
    month.lag0  = month.lag0,
    month.lag1  = month.lag1,
    year.lag0   = year.lag0,
    year.lag1   = year.lag1,
    wet.threshold = wet.threshold,
    extreme.threshold = extreme.threshold,
    month.list = month.list,
    sim.months = sim.months,
    sim.water.years = sim.water.years,
    year.idx = 1L,
    sim.start.year = 2000L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months)
  )

  idx <- which(sim.months == 7L)

  p00 <- res$p00_final[idx][1]
  p01 <- res$p01_final[idx][1]
  p02 <- res$p02_final[idx][1]

  testthat::expect_gt(p01, 0)
  testthat::expect_gt(p02, 0)

  testthat::expect_true(all(c(p00, p01, p02) >= 0))
  testthat::expect_true(all(c(p00, p01, p02) <= 1))

  testthat::expect_equal(p00 + p01 + p02, 1, tolerance = 1e-12)
})

testthat::test_that("monthly_markov_probs smoothing strength affects probabilities", {

  precip.lag1 <- c(rep(0, 10), rep(5, 2))  # mostly dry, a few wet
  precip.lag0 <- c(rep(0, 10), rep(5, 2))

  month.lag1 <- rep(6L, 12)
  month.lag0 <- rep(6L, 12)
  year.lag1  <- rep(2001L, 12)
  year.lag0  <- rep(2001L, 12)

  wet.threshold     <- rep(1, 12)
  extreme.threshold <- rep(10, 12)
  month.list <- 1:12

  sim.months <- rep(6L, 30)
  sim.water.years <- rep(2001L, 30)

  # Default smoothing (alpha default)
  res_low <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L, sim.start.year = 2001L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months)
  )

  idx <- which(sim.months == 6L)[1]
  p01_low <- res_low$p01_final[idx]

  # Increase sample size (reduces effective alpha_m = alpha/sqrt(N_m))
  precip.lag1b <- rep(precip.lag1, 5)
  precip.lag0b <- rep(precip.lag0, 5)
  month.lag1b  <- rep(month.lag1, 5)
  month.lag0b  <- rep(month.lag0, 5)
  year.lag1b   <- rep(year.lag1, 5)
  year.lag0b   <- rep(year.lag0, 5)

  res_high <- monthly_markov_probs(
    precip.lag0b, precip.lag1b,
    month.lag0b, month.lag1b,
    year.lag0b, year.lag1b,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L, sim.start.year = 2001L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months)
  )

  p01_high <- res_high$p01_final[idx]

  # With more data, smoothing influence should be weaker (closer to empirical 0)
  testthat::expect_gt(abs(p01_low - 0), abs(p01_high - 0))
})

testthat::test_that("monthly_markov_probs smoothing minimally affects well-sampled months", {

  set.seed(1)

  n <- 1000
  precip.lag1 <- rgamma(n, shape = 2, scale = 2)
  precip.lag0 <- rgamma(n, shape = 2, scale = 2)

  month.lag1 <- rep(8L, n)
  month.lag0 <- rep(8L, n)
  year.lag1  <- rep(1999L, n)
  year.lag0  <- rep(1999L, n)

  wet.threshold     <- rep(2, 12)
  extreme.threshold <- rep(6, 12)
  month.list <- 1:12

  sim.months <- rep(8L, 31)
  sim.water.years <- rep(1999L, 31)

  res <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L, sim.start.year = 1999L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months)
  )

  idx <- which(sim.months == 8L)[1]

  probs <- c(
    res$p00_final[idx], res$p01_final[idx], res$p02_final[idx],
    res$p10_final[idx], res$p11_final[idx], res$p12_final[idx]
  )

  testthat::expect_true(all(probs > 0.01))
  testthat::expect_true(all(probs < 0.99))
})

testthat::test_that("monthly_markov_probs produces no absorbing states", {

  precip.lag1 <- rep(0, 15)
  precip.lag0 <- rep(0, 15)

  month.lag1 <- rep(1L, 15)
  month.lag0 <- rep(1L, 15)
  year.lag1  <- rep(1980L, 15)
  year.lag0  <- rep(1980L, 15)

  wet.threshold     <- rep(1, 12)
  extreme.threshold <- rep(10, 12)
  month.list <- 1:12

  sim.months <- rep(1L, 31)
  sim.water.years <- rep(1980L, 31)

  res <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L, sim.start.year = 1980L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months)
  )

  idx <- which(sim.months == 1L)[1]

  testthat::expect_gt(sum(c(res$p00_final[idx], res$p01_final[idx], res$p02_final[idx]) > 0), 1)
  testthat::expect_gt(sum(c(res$p10_final[idx], res$p11_final[idx], res$p12_final[idx]) > 0), 1)
  testthat::expect_gt(sum(c(res$p20_final[idx], res$p21_final[idx], res$p22_final[idx]) > 0), 1)
})

testthat::test_that("monthly_markov_probs respects alpha argument", {

  precip.lag1 <- rep(0, 10)
  precip.lag0 <- rep(0, 10)
  month.lag1 <- rep(5L, 10)
  month.lag0 <- rep(5L, 10)
  year.lag1 <- rep(2000L, 10)
  year.lag0 <- rep(2000L, 10)

  wet.threshold     <- rep(1, 12)
  extreme.threshold <- rep(10, 12)
  month.list <- 1:12

  sim.months <- rep(5L, 31)
  sim.water.years <- rep(2000L, 31)

  res_nosmooth <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L, sim.start.year = 2000L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months),
    alpha = 0
  )

  res_smooth <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L, sim.start.year = 2000L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months),
    alpha = 0.5
  )

  idx <- which(sim.months == 5L)[1]

  testthat::expect_true(res_nosmooth$p01_final[idx] == 0)
  testthat::expect_gt(res_smooth$p01_final[idx], 0)
})

testthat::test_that("monthly_markov_probs alpha(N) decreases with sample size", {

  # Month 6: sparse
  precip.lag1_6 <- rep(0, 4)
  precip.lag0_6 <- rep(0, 4)

  # Month 7: dense
  precip.lag1_7 <- rep(0, 100)
  precip.lag0_7 <- rep(0, 100)

  precip.lag1 <- c(precip.lag1_6, precip.lag1_7)
  precip.lag0 <- c(precip.lag0_6, precip.lag0_7)

  month.lag1 <- c(rep(6L, 4), rep(7L, 100))
  month.lag0 <- month.lag1

  year.lag1 <- rep(2000L, length(precip.lag1))
  year.lag0 <- year.lag1

  wet.threshold     <- rep(1, 12)
  extreme.threshold <- rep(10, 12)
  month.list <- 1:12

  sim.months <- c(rep(6L, 30), rep(7L, 30))
  sim.water.years <- rep(2000L, length(sim.months))

  res <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L,
    sim.start.year = 2000L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months),
    alpha = 1
  )

  idx6 <- which(sim.months == 6L)[1]
  idx7 <- which(sim.months == 7L)[1]

  p01_6 <- res$p01_final[idx6]
  p01_7 <- res$p01_final[idx7]

  testthat::expect_gt(p01_6, p01_7)
})

testthat::test_that("monthly_markov_probs alpha(N) vanishes for large sample sizes", {

  set.seed(1)

  n <- 1000

  precip.lag1 <- c(rep(0, n / 2), rep(5, n / 2))
  precip.lag0 <- precip.lag1

  month.lag1 <- rep(9L, n)
  month.lag0 <- rep(9L, n)

  year.lag1 <- rep(1995L, n)
  year.lag0 <- year.lag1

  wet.threshold     <- rep(1, 12)
  extreme.threshold <- rep(10, 12)
  month.list <- 1:12

  sim.months <- rep(9L, 30)
  sim.water.years <- rep(1995L, 30)

  res <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L,
    sim.start.year = 1995L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months),
    alpha = 1
  )

  idx <- which(sim.months == 9L)[1]

  testthat::expect_gt(res$p00_final[idx], 0.95)
  testthat::expect_lt(res$p01_final[idx] + res$p02_final[idx], 0.05)
})

testthat::test_that("monthly_markov_probs rows sum to 1 for each origin state", {

  set.seed(42)

  n <- 300

  precip.lag1 <- rgamma(n, shape = 2, scale = 2)
  precip.lag0 <- rgamma(n, shape = 2, scale = 2)

  month.lag1 <- rep(3L, n)
  month.lag0 <- rep(3L, n)

  year.lag1 <- rep(2005L, n)
  year.lag0 <- rep(2005L, n)

  wet.threshold     <- rep(2, 12)
  extreme.threshold <- rep(8, 12)
  month.list <- 1:12

  sim.months <- rep(3L, 31)
  sim.water.years <- rep(2005L, 31)

  res <- monthly_markov_probs(
    precip.lag0, precip.lag1,
    month.lag0, month.lag1,
    year.lag0, year.lag1,
    wet.threshold, extreme.threshold,
    month.list,
    sim.months, sim.water.years,
    year.idx = 1L,
    sim.start.year = 2005L,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    sim.length = length(sim.months),
    alpha = 1
  )

  # pick any simulated index for March
  i <- which(sim.months == 3L)[1]

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

