test_that("valid index is returned unchanged", {
  prcp <- c(0, 5, 10, 20)

  idx <- get_result_index(2L, prcp)

  expect_identical(idx, 2L)
})


test_that("out-of-bounds index falls back to random valid index", {
  prcp <- c(0, 5, 10, 20)

  set.seed(123)
  idx <- get_result_index(10L, prcp)

  expect_true(idx %in% seq_along(prcp))
  expect_true(is.integer(idx) || is.numeric(idx))
})


test_that("negative or zero index falls back to random valid index", {
  prcp <- c(1, 2, 3)

  set.seed(42)
  idx0 <- get_result_index(0L, prcp)
  idxn <- get_result_index(-1L, prcp)

  expect_true(idx0 %in% seq_along(prcp))
  expect_true(idxn %in% seq_along(prcp))
})


test_that("NA index falls back to random valid index", {
  prcp <- c(3, 6, 9)

  set.seed(99)
  idx <- get_result_index(NA_integer_, prcp)

  expect_true(idx %in% seq_along(prcp))
})


test_that("length-zero RESULT falls back to random valid index", {
  prcp <- c(2, 4, 6)

  set.seed(7)
  idx <- get_result_index(integer(0), prcp)

  expect_true(idx %in% seq_along(prcp))
})


test_that("empty PRCP_TOMORROW always returns NA_integer_", {
  prcp_empty <- numeric(0)

  idx1 <- get_result_index(1L, prcp_empty)
  idx2 <- get_result_index(NA_integer_, prcp_empty)
  idx3 <- get_result_index(100L, prcp_empty)

  expect_identical(idx1, NA_integer_)
  expect_identical(idx2, NA_integer_)
  expect_identical(idx3, NA_integer_)
})


test_that("returned index always lies within bounds when PRCP_TOMORROW is non-empty", {
  prcp <- rnorm(10)

  set.seed(2024)
  idxs <- replicate(
    100,
    get_result_index(sample(c(NA, -5, 0, 15), 1), prcp)
  )

  expect_true(all(idxs %in% seq_along(prcp)))
})
