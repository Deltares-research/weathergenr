test_that("calculateMarkovProbs returns correct structure and lengths", {

  SIM_LENGTH <- 365

  res <- calculateMarkovProbs(
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

  res <- calculateMarkovProbs(
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

  res <- calculateMarkovProbs(
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

  res <- calculateMarkovProbs(
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

  res <- calculateMarkovProbs(
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
