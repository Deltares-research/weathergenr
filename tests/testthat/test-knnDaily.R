library(testthat)

test_that("knnDaily returns an integer index within correct bounds", {
  set.seed(123)
  cur_sim_PRCP <- 7
  cur_sim_TEMP <- 23
  PRCP_TODAY <- seq(0, 20, length.out = 50)
  TEMP_TODAY <- seq(15, 30, length.out = 50)
  k <- 5
  sd_monthly_PRCP <- 5
  sd_monthly_TEMP <- 2
  mean_monthly_PRCP <- 10
  mean_monthly_TEMP <- 20
  k1 <- 1
  count <- 1
  seed <- 100

  result <- knnDaily(
    cur_sim_PRCP = cur_sim_PRCP,
    cur_sim_TEMP = cur_sim_TEMP,
    PRCP_TODAY = PRCP_TODAY,
    TEMP_TODAY = TEMP_TODAY,
    k = k,
    sd_monthly_PRCP = sd_monthly_PRCP,
    sd_monthly_TEMP = sd_monthly_TEMP,
    mean_monthly_PRCP = mean_monthly_PRCP,
    mean_monthly_TEMP = mean_monthly_TEMP,
    k1 = k1,
    count = count,
    seed = seed
  )
  expect_type(result, "integer")
  expect_true(result >= 1 && result <= length(PRCP_TODAY))
})

test_that("knnDaily is reproducible with the same seed and settings", {
  cur_sim_PRCP <- 9
  cur_sim_TEMP <- 20
  PRCP_TODAY <- runif(80, 0, 30)
  TEMP_TODAY <- runif(80, 10, 35)
  k <- 7
  sd_monthly_PRCP <- 6
  sd_monthly_TEMP <- 4
  mean_monthly_PRCP <- 15
  mean_monthly_TEMP <- 22
  k1 <- 2
  count <- 4
  seed <- 77

  res1 <- knnDaily(
    cur_sim_PRCP, cur_sim_TEMP,
    PRCP_TODAY, TEMP_TODAY,
    k, sd_monthly_PRCP, sd_monthly_TEMP,
    mean_monthly_PRCP, mean_monthly_TEMP,
    k1, count, seed
  )
  res2 <- knnDaily(
    cur_sim_PRCP, cur_sim_TEMP,
    PRCP_TODAY, TEMP_TODAY,
    k, sd_monthly_PRCP, sd_monthly_TEMP,
    mean_monthly_PRCP, mean_monthly_TEMP,
    k1, count, seed
  )
  expect_identical(res1, res2)
})

test_that("knnDaily works when k = 1", {
  cur_sim_PRCP <- 10
  cur_sim_TEMP <- 18
  PRCP_TODAY <- runif(12, 0, 15)
  TEMP_TODAY <- runif(12, 10, 25)
  k <- 1
  sd_monthly_PRCP <- 2
  sd_monthly_TEMP <- 1
  mean_monthly_PRCP <- 5
  mean_monthly_TEMP <- 14
  k1 <- 3
  count <- 2
  seed <- 888

  result <- knnDaily(
    cur_sim_PRCP, cur_sim_TEMP,
    PRCP_TODAY, TEMP_TODAY,
    k, sd_monthly_PRCP, sd_monthly_TEMP,
    mean_monthly_PRCP, mean_monthly_TEMP,
    k1, count, seed
  )
  expect_type(result, "integer")
  expect_true(result >= 1 && result <= length(PRCP_TODAY))
})
