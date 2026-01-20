

test_that("valid index is returned unchanged", {
  precip <- c(0, 5, 10, 20)

  idx <- get_result_index(2L, precip)

  expect_identical(idx, 2L)
})


test_that("out-of-bounds index falls back to random valid index", {
  precip <- c(0, 5, 10, 20)

  set.seed(123)
  idx <- get_result_index(10L, precip)

  expect_true(idx %in% seq_along(precip))
  expect_true(is.integer(idx) || is.numeric(idx))
})


test_that("negative or zero index falls back to random valid index", {
  precip <- c(1, 2, 3)

  set.seed(42)
  idx0 <- get_result_index(0L, precip)
  idxn <- get_result_index(-1L, precip)

  expect_true(idx0 %in% seq_along(precip))
  expect_true(idxn %in% seq_along(precip))
})


test_that("NA index falls back to random valid index", {
  precip <- c(3, 6, 9)

  set.seed(99)
  idx <- get_result_index(NA_integer_, precip)

  expect_true(idx %in% seq_along(precip))
})


test_that("length-zero RESULT falls back to random valid index", {
  precip <- c(2, 4, 6)

  set.seed(7)
  idx <- get_result_index(integer(0), precip)

  expect_true(idx %in% seq_along(precip))
})


test_that("empty precip_TOMORROW always returns NA_integer_", {
  precip_empty <- numeric(0)

  idx1 <- get_result_index(1L, precip_empty)
  idx2 <- get_result_index(NA_integer_, precip_empty)
  idx3 <- get_result_index(100L, precip_empty)

  expect_identical(idx1, NA_integer_)
  expect_identical(idx2, NA_integer_)
  expect_identical(idx3, NA_integer_)
})


test_that("returned index always lies within bounds when precip_TOMORROW is non-empty", {
  precip <- rnorm(10)

  set.seed(2024)
  idxs <- replicate(
    100,
    get_result_index(sample(c(NA, -5, 0, 15), 1), precip)
  )

  expect_true(all(idxs %in% seq_along(precip)))
})
