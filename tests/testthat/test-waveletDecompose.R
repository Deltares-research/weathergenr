library(testthat)

# Example synthetic data
set.seed(42)
x <- as.numeric(scale(rnorm(40))) # Length-40 "annual" series

# Assume significant periods 2 and 4
sig_periods <- list(c(2), c(4))

test_that("waveletDecompose returns tibble with expected columns", {
  result <- waveletDecompose(variable = x, signif.periods = sig_periods, plot = FALSE)

  expect_true(is.data.frame(result) || inherits(result, "tbl"))
  expect_true("NOISE" %in% names(result))
  expect_true(any(grepl("^Component_", names(result))))
})

test_that("Components and NOISE have correct length", {
  result <- waveletDecompose(variable = x, signif.periods = sig_periods, plot = FALSE)
  n <- length(x)
  expect_true(all(vapply(result, length, integer(1)) == n))
})

test_that("waveletDecompose errors for NULL variable", {
  expect_error(
    waveletDecompose(variable = NULL, signif.periods = sig_periods, plot = FALSE),
    "NULL|missing|variable"
  )
})

# test_that("waveletDecompose errors for NULL signif.periods", {
#   x <- rnorm(100)
#   expect_error(waveletDecompose(variable = x, signif.periods = NULL, plot = FALSE), "must not be NULL")
# })

# test_that("waveletDecompose returns all zeros if no significant periods", {
#   result <- waveletDecompose(variable = rep(0, 40), signif.periods = list(), plot = FALSE)
#   expect_true(all(result$NOISE == 0))
# })


# Edge case: length of variable and signif.periods mismatch
test_that("waveletDecompose returns correct structure even with 1 period", {
  x <- as.numeric(scale(rnorm(10)))
  sig_periods <- list(2)
  result <- waveletDecompose(variable = x, signif.periods = sig_periods, plot = FALSE)
  expect_true("Component_1" %in% names(result))
})
