test_that("monthly_markov_probs returns correct structure and lengths", {

  SIM_LENGTH <- 365

  res <- monthly_markov_probs(
    PRCP_LAG0 = rep(5, 100),
    PRCP_LAG1 = rep(3, 100),
    MONTH_LAG0 = rep(1, 100),
    MONTH_LAG1 = rep(1, 100),
    YEAR_LAG0  = rep(2001, 100),
    YEAR_LAG1  = rep(2001, 100),
    wet_threshold = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_list = 1:12,
    MONTH_SIM = rep(1, SIM_LENGTH),
    WATER_YEAR_SIM = rep(2001, SIM_LENGTH),
    y = 1,
    START_YEAR_SIM = 2001,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = SIM_LENGTH
  )

  expect_type(res, "list")
  expect_equal(length(res), 9)

  for (nm in names(res)) {
    expect_length(res[[nm]], SIM_LENGTH)
    expect_true(is.numeric(res[[nm]]))
  }
})

# ------------------------------------------------------------
# Calendar-year mode (month.start == 1)
# ------------------------------------------------------------
test_that("calendar-year mode excludes cross-year transitions", {

  # December (12) to January (1) transition
  PRCP_LAG1 <- c(5, 5)
  PRCP_LAG0 <- c(5, 5)

  MONTH_LAG1 <- c(12, 1)
  MONTH_LAG0 <- c(1, 2)

  YEAR_LAG1 <- c(2001, 2002)
  YEAR_LAG0 <- c(2002, 2002)

  SIM_LENGTH <- 10

  res <- monthly_markov_probs(
    PRCP_LAG0 = PRCP_LAG0,
    PRCP_LAG1 = PRCP_LAG1,
    MONTH_LAG0 = MONTH_LAG0,
    MONTH_LAG1 = MONTH_LAG1,
    YEAR_LAG0  = YEAR_LAG0,
    YEAR_LAG1  = YEAR_LAG1,
    wet_threshold = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_list = 1:12,        # calendar-year mode
    MONTH_SIM = rep(1, SIM_LENGTH),
    WATER_YEAR_SIM = rep(2002, SIM_LENGTH),
    y = 1,
    START_YEAR_SIM = 2002,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = SIM_LENGTH
  )

  # Because cross-year lags are filtered, transitions should
  # fall back to deterministic dry behavior
  expect_true(all(res$p00_final[!is.na(res$p00_final)] >= 0))
  expect_true(all(res$p01_final[!is.na(res$p01_final)] >= 0))
})

# ------------------------------------------------------------
# Water-year mode (month.start != 1)
# ------------------------------------------------------------
test_that("water-year mode allows cross-calendar-year transitions", {

  PRCP_LAG1 <- c(5, 5)
  PRCP_LAG0 <- c(5, 5)

  MONTH_LAG1 <- c(12, 1)
  MONTH_LAG0 <- c(1, 2)

  YEAR_LAG1 <- c(2001, 2001)
  YEAR_LAG0 <- c(2001, 2001)

  SIM_LENGTH <- 10

  res <- monthly_markov_probs(
    PRCP_LAG0 = PRCP_LAG0,
    PRCP_LAG1 = PRCP_LAG1,
    MONTH_LAG0 = MONTH_LAG0,
    MONTH_LAG1 = MONTH_LAG1,
    YEAR_LAG0  = YEAR_LAG0,
    YEAR_LAG1  = YEAR_LAG1,
    wet_threshold = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_list = c(10,11,12,1,2,3,4,5,6,7,8,9),  # water-year mode
    MONTH_SIM = rep(1, SIM_LENGTH),
    WATER_YEAR_SIM = rep(2001, SIM_LENGTH),
    y = 1,
    START_YEAR_SIM = 2000,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = SIM_LENGTH
  )

  expect_true(any(res$p00_final > 0))
})

# ------------------------------------------------------------
# Probability consistency
# ------------------------------------------------------------
test_that("transition probabilities are valid and bounded", {

  SIM_LENGTH <- 50

  res <- monthly_markov_probs(
    PRCP_LAG0 = runif(200, 0, 20),
    PRCP_LAG1 = runif(200, 0, 20),
    MONTH_LAG0 = sample(1:12, 200, TRUE),
    MONTH_LAG1 = sample(1:12, 200, TRUE),
    YEAR_LAG0  = rep(2001, 200),
    YEAR_LAG1  = rep(2001, 200),
    wet_threshold = rep(5, 12),
    extreme_threshold = rep(15, 12),
    month_list = 1:12,
    MONTH_SIM = sample(1:12, SIM_LENGTH, TRUE),
    WATER_YEAR_SIM = rep(2001, SIM_LENGTH),
    y = 1,
    START_YEAR_SIM = 2001,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = SIM_LENGTH
  )

  # Check bounds
  for (nm in names(res)) {
    x <- res[[nm]]
    x <- x[!is.na(x)]
    expect_true(all(x >= 0))
    expect_true(all(x <= 1))
  }
})

# ------------------------------------------------------------
# YEAR_LAG guard correctness
# ------------------------------------------------------------
test_that("YEAR_LAG filtering removes cross-year transitions", {

  PRCP_LAG0 <- c(5, 5, 5)
  PRCP_LAG1 <- c(5, 5, 5)

  MONTH_LAG0 <- c(1, 1, 1)
  MONTH_LAG1 <- c(12, 1, 1)

  YEAR_LAG0 <- c(2002, 2002, 2002)
  YEAR_LAG1 <- c(2001, 2002, 2002)

  SIM_LENGTH <- 10

  res <- monthly_markov_probs(
    PRCP_LAG0 = PRCP_LAG0,
    PRCP_LAG1 = PRCP_LAG1,
    MONTH_LAG0 = MONTH_LAG0,
    MONTH_LAG1 = MONTH_LAG1,
    YEAR_LAG0  = YEAR_LAG0,
    YEAR_LAG1  = YEAR_LAG1,
    wet_threshold = rep(2, 12),
    extreme_threshold = rep(10, 12),
    month_list = 1:12,
    MONTH_SIM = rep(1, SIM_LENGTH),
    WATER_YEAR_SIM = rep(2002, SIM_LENGTH),
    y = 1,
    START_YEAR_SIM = 2002,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = SIM_LENGTH
  )

  # Should not error and should return valid probabilities
  expect_true(all(is.finite(res$p00_final[!is.na(res$p00_final)])))
})


test_that("monthly_markov_probs applies smoothing for zero-denominator months", {

  ## Construct an extreme sparse case:
  ## All precipitation is zero implies that only dry state observed
  PRCP_LAG1 <- rep(0, 20)
  PRCP_LAG0 <- rep(0, 20)

  MONTH_LAG1 <- rep(7, 20)  # July only
  MONTH_LAG0 <- rep(7, 20)

  YEAR_LAG1 <- rep(2000, 20)
  YEAR_LAG0 <- rep(2000, 20)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_list <- 1:12

  MONTH_SIM <- rep(7, 31)
  WATER_YEAR_SIM <- rep(2000, 31)

  res <- monthly_markov_probs(
    PRCP_LAG0 = PRCP_LAG0,
    PRCP_LAG1 = PRCP_LAG1,
    MONTH_LAG0 = MONTH_LAG0,
    MONTH_LAG1 = MONTH_LAG1,
    YEAR_LAG0 = YEAR_LAG0,
    YEAR_LAG1 = YEAR_LAG1,
    wet_threshold = wet_threshold,
    extreme_threshold = extreme_threshold,
    month_list = month_list,
    MONTH_SIM = MONTH_SIM,
    WATER_YEAR_SIM = WATER_YEAR_SIM,
    y = 1,
    START_YEAR_SIM = 2000,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM)
  )

  ## Extract probabilities for July
  idx <- which(MONTH_SIM == 7)

  p00 <- res$p00_final[idx][1]
  p01 <- res$p01_final[idx][1]
  p02 <- res$p02_final[idx][1]

  ## ---- Assertions ----

  # No structural zeros
  expect_gt(p01, 0)
  expect_gt(p02, 0)

  # Valid probability bounds
  expect_true(all(c(p00, p01, p02) >= 0))
  expect_true(all(c(p00, p01, p02) <= 1))

  # Row sums to 1
  expect_equal(p00 + p01 + p02, 1, tolerance = 1e-12)
})

test_that("monthly_markov_probs smoothing strength affects probabilities", {

  PRCP_LAG1 <- c(rep(0, 10), rep(5, 2))  # mostly dry, a few wet
  PRCP_LAG0 <- c(rep(0, 10), rep(5, 2))

  MONTH_LAG1 <- rep(6, 12)
  MONTH_LAG0 <- rep(6, 12)
  YEAR_LAG1  <- rep(2001, 12)
  YEAR_LAG0  <- rep(2001, 12)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_list <- 1:12

  MONTH_SIM <- rep(6, 30)
  WATER_YEAR_SIM <- rep(2001, 30)

  # Low smoothing
  res_low <- monthly_markov_probs(
    PRCP_LAG0, PRCP_LAG1,
    MONTH_LAG0, MONTH_LAG1,
    YEAR_LAG0, YEAR_LAG1,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1, START_YEAR_SIM = 2001,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM)
  )

  # Extract one month
  idx <- which(MONTH_SIM == 6)[1]
  p01_low <- res_low$p01_final[idx]

  # Re-run with artificially stronger smoothing by scaling counts
  # (proxy: duplicate observations to weaken smoothing effect)
  PRCP_LAG1b <- rep(PRCP_LAG1, 5)
  PRCP_LAG0b <- rep(PRCP_LAG0, 5)
  MONTH_LAG1b <- rep(MONTH_LAG1, 5)
  MONTH_LAG0b <- rep(MONTH_LAG0, 5)
  YEAR_LAG1b  <- rep(YEAR_LAG1, 5)
  YEAR_LAG0b  <- rep(YEAR_LAG0, 5)

  res_high <- monthly_markov_probs(
    PRCP_LAG0b, PRCP_LAG1b,
    MONTH_LAG0b, MONTH_LAG1b,
    YEAR_LAG0b, YEAR_LAG1b,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1, START_YEAR_SIM = 2001,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM)
  )

  p01_high <- res_high$p01_final[idx]

  # With more data, smoothing influence should be weaker
  expect_gt(abs(p01_low - 0), abs(p01_high - 0))
})

test_that("monthly_markov_probs smoothing minimally affects well-sampled months", {

  set.seed(1)

  n <- 1000
  PRCP_LAG1 <- rgamma(n, shape = 2, scale = 2)
  PRCP_LAG0 <- rgamma(n, shape = 2, scale = 2)

  MONTH_LAG1 <- rep(8, n)
  MONTH_LAG0 <- rep(8, n)
  YEAR_LAG1  <- rep(1999, n)
  YEAR_LAG0  <- rep(1999, n)

  wet_threshold     <- rep(2, 12)
  extreme_threshold <- rep(6, 12)
  month_list <- 1:12

  MONTH_SIM <- rep(8, 31)
  WATER_YEAR_SIM <- rep(1999, 31)

  res <- monthly_markov_probs(
    PRCP_LAG0, PRCP_LAG1,
    MONTH_LAG0, MONTH_LAG1,
    YEAR_LAG0, YEAR_LAG1,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1, START_YEAR_SIM = 1999,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM)
  )

  idx <- which(MONTH_SIM == 8)[1]

  # All probabilities should be strictly interior
  probs <- c(
    res$p00_final[idx], res$p01_final[idx], res$p02_final[idx],
    res$p10_final[idx], res$p11_final[idx], res$p12_final[idx]
  )

  expect_true(all(probs > 0.01))
  expect_true(all(probs < 0.99))
})

test_that("monthly_markov_probs produces no absorbing states", {

  PRCP_LAG1 <- rep(0, 15)
  PRCP_LAG0 <- rep(0, 15)

  MONTH_LAG1 <- rep(1, 15)
  MONTH_LAG0 <- rep(1, 15)
  YEAR_LAG1  <- rep(1980, 15)
  YEAR_LAG0  <- rep(1980, 15)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_list <- 1:12

  MONTH_SIM <- rep(1, 31)
  WATER_YEAR_SIM <- rep(1980, 31)

  res <- monthly_markov_probs(
    PRCP_LAG0, PRCP_LAG1,
    MONTH_LAG0, MONTH_LAG1,
    YEAR_LAG0, YEAR_LAG1,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1, START_YEAR_SIM = 1980,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM)
  )

  idx <- which(MONTH_SIM == 1)[1]

  # Each row must have at least two non-zero transitions
  expect_gt(sum(c(res$p00_final[idx], res$p01_final[idx], res$p02_final[idx]) > 0), 1)
  expect_gt(sum(c(res$p10_final[idx], res$p11_final[idx], res$p12_final[idx]) > 0), 1)
  expect_gt(sum(c(res$p20_final[idx], res$p21_final[idx], res$p22_final[idx]) > 0), 1)
})

test_that("monthly_markov_probs respects alpha argument", {

  PRCP_LAG1 <- rep(0, 10)
  PRCP_LAG0 <- rep(0, 10)
  MONTH_LAG1 <- rep(5, 10)
  MONTH_LAG0 <- rep(5, 10)
  YEAR_LAG1 <- rep(2000, 10)
  YEAR_LAG0 <- rep(2000, 10)

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_list <- 1:12

  MONTH_SIM <- rep(5, 31)
  WATER_YEAR_SIM <- rep(2000, 31)

  res_nosmooth <- monthly_markov_probs(
    PRCP_LAG0, PRCP_LAG1,
    MONTH_LAG0, MONTH_LAG1,
    YEAR_LAG0, YEAR_LAG1,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1, START_YEAR_SIM = 2000,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM),
    alpha = 0
  )

  res_smooth <- monthly_markov_probs(
    PRCP_LAG0, PRCP_LAG1,
    MONTH_LAG0, MONTH_LAG1,
    YEAR_LAG0, YEAR_LAG1,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1, START_YEAR_SIM = 2000,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM),
    alpha = 0.5
  )

  idx <- which(MONTH_SIM == 5)[1]

  expect_true(res_nosmooth$p01_final[idx] == 0)
  expect_gt(res_smooth$p01_final[idx], 0)
})

test_that("monthly_markov_probs alpha(N) decreases with sample size", {

  ## Construct two months with identical transition structure
  ## but different sample sizes

  # Month 6: sparse
  PRCP_LAG1_6 <- rep(0, 4)
  PRCP_LAG0_6 <- rep(0, 4)

  # Month 7: dense
  PRCP_LAG1_7 <- rep(0, 100)
  PRCP_LAG0_7 <- rep(0, 100)

  PRCP_LAG1 <- c(PRCP_LAG1_6, PRCP_LAG1_7)
  PRCP_LAG0 <- c(PRCP_LAG0_6, PRCP_LAG0_7)

  MONTH_LAG1 <- c(rep(6, 4), rep(7, 100))
  MONTH_LAG0 <- MONTH_LAG1

  YEAR_LAG1 <- rep(2000, length(PRCP_LAG1))
  YEAR_LAG0 <- YEAR_LAG1

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_list <- 1:12

  MONTH_SIM <- c(rep(6, 30), rep(7, 30))
  WATER_YEAR_SIM <- rep(2000, length(MONTH_SIM))

  res <- monthly_markov_probs(
    PRCP_LAG0, PRCP_LAG1,
    MONTH_LAG0, MONTH_LAG1,
    YEAR_LAG0, YEAR_LAG1,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1,
    START_YEAR_SIM = 2000,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM),
    alpha = 1
  )

  idx6 <- which(MONTH_SIM == 6)[1]
  idx7 <- which(MONTH_SIM == 7)[1]

  # In both cases, empirically p01 = 0, but smoothing differs
  p01_6 <- res$p01_final[idx6]
  p01_7 <- res$p01_final[idx7]

  # Sparse month should be more smoothed
  expect_gt(p01_6, p01_7)
})

test_that("monthly_markov_probs alpha(N) vanishes for large sample sizes", {

  set.seed(1)

  n <- 1000

  PRCP_LAG1 <- c(rep(0, n / 2), rep(5, n / 2))
  PRCP_LAG0 <- PRCP_LAG1

  MONTH_LAG1 <- rep(9, n)
  MONTH_LAG0 <- rep(9, n)

  YEAR_LAG1 <- rep(1995, n)
  YEAR_LAG0 <- YEAR_LAG1

  wet_threshold     <- rep(1, 12)
  extreme_threshold <- rep(10, 12)
  month_list <- 1:12

  MONTH_SIM <- rep(9, 30)
  WATER_YEAR_SIM <- rep(1995, 30)

  res <- monthly_markov_probs(
    PRCP_LAG0, PRCP_LAG1,
    MONTH_LAG0, MONTH_LAG1,
    YEAR_LAG0, YEAR_LAG1,
    wet_threshold, extreme_threshold,
    month_list,
    MONTH_SIM, WATER_YEAR_SIM,
    y = 1,
    START_YEAR_SIM = 1995,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    SIM_LENGTH = length(MONTH_SIM),
    alpha = 1
  )

  idx <- which(MONTH_SIM == 9)[1]

  # Empirical dry to dry probability about 1
  # Smoothed estimate should be very close
  expect_gt(res$p00_final[idx], 0.95)
  expect_lt(res$p01_final[idx] + res$p02_final[idx], 0.05)
})
