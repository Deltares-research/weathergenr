library(testthat)

library(testthat)

test_that("getPI returns valid stationary distribution for a standard Markov transition matrix", {
  result <- getPI(
    p00 = 0.6, p01 = 0.3, p02 = 0.1,
    p10 = 0.2, p11 = 0.7, p12 = 0.1,
    p20 = 0.3, p21 = 0.2, p22 = 0.5
  )
  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("getPI returns valid stationary distribution for high persistence", {
  result <- getPI(
    p00 = 0.9, p01 = 0.05, p02 = 0.05,
    p10 = 0.1, p11 = 0.85, p12 = 0.05,
    p20 = 0.2, p21 = 0.2, p22 = 0.6
  )
  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("getPI works with nearly deterministic chains", {
  result <- getPI(
    p00 = 0.99, p01 = 0.01, p02 = 0,
    p10 = 0, p11 = 0.99, p12 = 0.01,
    p20 = 0.01, p21 = 0, p22 = 0.99
  )
  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})
