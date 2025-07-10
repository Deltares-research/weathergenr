library(testthat)

# Minimal reproducible input generator for ARIMA
make_components <- function(n_years = 20, n_comp = 2, seed = 123) {
  set.seed(seed)
  as.data.frame(matrix(rnorm(n_years * n_comp), ncol = n_comp))
}

test_that("waveletARIMA returns matrix with correct dimensions", {
  # Simulate two components, 10 years, 5 traces
  comp <- make_components(n_years = 10, n_comp = 2)
  out <- waveletARIMA(
    wavelet.components = comp,
    sim.year.num = 10,
    sim.num = 5,
    seed = 42
  )
  expect_true(is.matrix(out) || is.data.frame(out))
  expect_equal(dim(out), c(10, 5))
})

test_that("waveletARIMA output is reproducible with the same seed", {
  comp <- make_components(n_years = 12, n_comp = 2)
  res1 <- waveletARIMA(comp, sim.year.num = 12, sim.num = 3, seed = 99)
  res2 <- waveletARIMA(comp, sim.year.num = 12, sim.num = 3, seed = 99)
  expect_equal(res1, res2)
})

test_that("waveletARIMA output varies with different seeds", {
  comp <- make_components(n_years = 12, n_comp = 2)
  res1 <- waveletARIMA(comp, sim.year.num = 12, sim.num = 3, seed = 1)
  res2 <- waveletARIMA(comp, sim.year.num = 12, sim.num = 3, seed = 2)
  expect_false(isTRUE(all.equal(res1, res2)))
})

test_that("waveletARIMA works for a single component (column)", {
  comp <- make_components(n_years = 8, n_comp = 1)
  out <- waveletARIMA(comp, sim.year.num = 8, sim.num = 4, seed = 101)
  expect_true(is.matrix(out) || is.data.frame(out))
  expect_equal(dim(out), c(8, 4))
})


test_that("waveletARIMA returns error if input is NULL", {
  expect_error(
    waveletARIMA(NULL, sim.year.num = 10, sim.num = 2, seed = 10),
    regexp = "NULL|missing|Error"
  )
})

# Edge case: sim.num = 1
test_that("waveletARIMA works for sim.num = 1", {
  comp <- make_components(n_years = 7, n_comp = 2)
  out <- waveletARIMA(comp, sim.year.num = 7, sim.num = 1, seed = 123)
  expect_true(is.matrix(out) || is.data.frame(out))
  expect_equal(dim(out), c(7, 1))
})
