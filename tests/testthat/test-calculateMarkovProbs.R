test_that("calculateMarkovProbs returns correct output structure and lengths", {
  set.seed(42)
  n <- 100
  PRCP_LAG1 <- runif(n, 0, 20)
  PRCP_LAG0 <- runif(n, 0, 20)
  MONTH_LAG1 <- rep(1, n)
  MONTH_LAG0 <- rep(1, n)
  wet_threshold <- rep(5, 12)
  extreme_threshold <- rep(15, 12)
  month_list <- 1:12
  MONTH_SIM <- rep(1, n)
  WATER_YEAR_SIM <- rep(2020, n)
  y <- 0
  START_YEAR_SIM <- 2020
  dry.spell.change <- rep(1, 12)
  wet.spell.change <- rep(1, 12)
  SIM_LENGTH <- n

  result <- calculateMarkovProbs(
    PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
    wet_threshold, extreme_threshold,
    month_list, MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
    dry.spell.change, wet.spell.change, SIM_LENGTH
  )
  expect_type(result, "list")
  expect_length(result, 9)
  expect_true(all(vapply(result, length, integer(1)) == n))
})

test_that("Transition probabilities sum to 1 or less for each index", {
  # Generate reproducible dummy data for a single month
  n <- 100
  PRCP_LAG1 <- runif(n, 0, 20)
  PRCP_LAG0 <- runif(n, 0, 20)
  MONTH_LAG1 <- rep(1, n)
  MONTH_LAG0 <- rep(1, n)
  wet_threshold <- rep(5, 12)
  extreme_threshold <- rep(15, 12)
  month_list <- 1:12
  MONTH_SIM <- rep(1, n)
  WATER_YEAR_SIM <- rep(2020, n)
  y <- 0
  START_YEAR_SIM <- 2020
  dry.spell.change <- rep(1, 12)
  wet.spell.change <- rep(1, 12)
  SIM_LENGTH <- n

  result <- calculateMarkovProbs(
    PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
    wet_threshold, extreme_threshold,
    month_list, MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
    dry.spell.change, wet.spell.change, SIM_LENGTH
  )
  # For each time step, one of the p00, p01, p02 (dry state), etc, should sum to <=1 (could be <1 if denom is 0)
  for (i in 1:n) {
    expect_true(result$p00_final[i] + result$p01_final[i] + result$p02_final[i] <= 1 + 1e-8)
    expect_true(result$p10_final[i] + result$p11_final[i] + result$p12_final[i] <= 1 + 1e-8)
    expect_true(result$p20_final[i] + result$p21_final[i] + result$p22_final[i] <= 1 + 1e-8)
  }
})

test_that("Handles no matching month (returns zeros for that month)", {
  n <- 50
  # All months are set to February, but we only check for January
  PRCP_LAG1 <- runif(n, 0, 20)
  PRCP_LAG0 <- runif(n, 0, 20)
  MONTH_LAG1 <- rep(2, n)
  MONTH_LAG0 <- rep(2, n)
  wet_threshold <- rep(5, 12)
  extreme_threshold <- rep(15, 12)
  month_list <- 1:12
  MONTH_SIM <- rep(1, n) # Only January in simulation
  WATER_YEAR_SIM <- rep(2020, n)
  y <- 0
  START_YEAR_SIM <- 2020
  dry.spell.change <- rep(1, 12)
  wet.spell.change <- rep(1, 12)
  SIM_LENGTH <- n

  result <- calculateMarkovProbs(
    PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
    wet_threshold, extreme_threshold,
    month_list, MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
    dry.spell.change, wet.spell.change, SIM_LENGTH
  )
  # All probabilities should be zero since no month matches
  expect_true(all(unlist(result) == 0))
})

test_that("Works for degenerate states (all dry or all wet)", {
  n <- 30
  PRCP_LAG1 <- rep(0, n)
  PRCP_LAG0 <- rep(0, n)
  MONTH_LAG1 <- rep(1, n)
  MONTH_LAG0 <- rep(1, n)
  wet_threshold <- rep(5, 12)
  extreme_threshold <- rep(15, 12)
  month_list <- 1:12
  MONTH_SIM <- rep(1, n)
  WATER_YEAR_SIM <- rep(2020, n)
  y <- 0
  START_YEAR_SIM <- 2020
  dry.spell.change <- rep(1, 12)
  wet.spell.change <- rep(1, 12)
  SIM_LENGTH <- n

  result <- calculateMarkovProbs(
    PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
    wet_threshold, extreme_threshold,
    month_list, MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
    dry.spell.change, wet.spell.change, SIM_LENGTH
  )
  # At least some probability mass should be on p00_final (dry->dry)
  expect_true(any(result$p00_final > 0))
})

test_that("Handles short time series (n = 1)", {
  n <- 1
  PRCP_LAG1 <- 1
  PRCP_LAG0 <- 1
  MONTH_LAG1 <- 1
  MONTH_LAG0 <- 1
  wet_threshold <- rep(5, 12)
  extreme_threshold <- rep(15, 12)
  month_list <- 1:12
  MONTH_SIM <- 1
  WATER_YEAR_SIM <- 2020
  y <- 0
  START_YEAR_SIM <- 2020
  dry.spell.change <- rep(1, 12)
  wet.spell.change <- rep(1, 12)
  SIM_LENGTH <- 1

  result <- calculateMarkovProbs(
    PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
    wet_threshold, extreme_threshold,
    month_list, MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
    dry.spell.change, wet.spell.change, SIM_LENGTH
  )
  expect_length(result$p00_final, 1)
  expect_length(result$p11_final, 1)
})
