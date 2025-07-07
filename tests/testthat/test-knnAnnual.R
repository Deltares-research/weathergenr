library(testthat)

test_that("knnAnnual returns the correct length and type", {
  set.seed(42)
  sim_annual_prcp <- 850
  ANNUAL_PRCP <- c(800, 900, 750, 860, 880, 840, 830, 870, 810, 890)
  WATER_YEAR_A <- 2001:2010
  kk <- 5
  k1 <- 2
  y <- 3
  y_sample_size <- 7
  seed <- 123

  result <- knnAnnual(sim_annual_prcp, ANNUAL_PRCP, WATER_YEAR_A,
                      kk, k1, y, y_sample_size, seed)

  expect_length(result, y_sample_size)
  expect_type(result, "integer")
  expect_true(all(result %in% WATER_YEAR_A))
})

test_that("knnAnnual is reproducible given a seed", {
  sim_annual_prcp <- 820
  ANNUAL_PRCP <- c(810, 820, 800, 830, 815)
  WATER_YEAR_A <- 2011:2015
  kk <- 3
  k1 <- 5
  y <- 2
  y_sample_size <- 6
  seed <- 42

  result1 <- knnAnnual(sim_annual_prcp, ANNUAL_PRCP, WATER_YEAR_A,
                       kk, k1, y, y_sample_size, seed)
  result2 <- knnAnnual(sim_annual_prcp, ANNUAL_PRCP, WATER_YEAR_A,
                       kk, k1, y, y_sample_size, seed)
  expect_identical(result1, result2)
})

test_that("knnAnnual works with kk = length(ANNUAL_PRCP)", {
  sim_annual_prcp <- 855
  ANNUAL_PRCP <- c(850, 860, 870, 840, 855)
  WATER_YEAR_A <- 1996:2000
  kk <- length(ANNUAL_PRCP)
  k1 <- 1
  y <- 1
  y_sample_size <- 10
  result <- knnAnnual(sim_annual_prcp, ANNUAL_PRCP, WATER_YEAR_A,
                      kk, k1, y, y_sample_size)
  expect_length(result, y_sample_size)
  expect_true(all(result %in% WATER_YEAR_A))
})

test_that("knnAnnual returns values even if kk = 1", {
  sim_annual_prcp <- 900
  ANNUAL_PRCP <- c(850, 900, 870)
  WATER_YEAR_A <- 2010:2012
  kk <- 1
  k1 <- 1
  y <- 1
  y_sample_size <- 4
  result <- knnAnnual(sim_annual_prcp, ANNUAL_PRCP, WATER_YEAR_A,
                      kk, k1, y, y_sample_size)
  expect_length(result, y_sample_size)
  expect_true(all(result %in% WATER_YEAR_A))
})
