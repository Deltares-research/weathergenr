testthat::test_that("get_state_indices returns correct transitions for simple example", {

  prcp <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
  candidate.idx <- 1:(length(prcp) - 1)
  wet.thr <- 1
  extreme.thr <- 20

  # dry (0) -> wet (1)
  # transitions at:
  # idx 2: 0 -> 5
  # idx 6: 0 -> 2
  res01 <- get_state_indices(0, 1, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res01, c(2L, 6L))

  # wet (1) -> extreme (2)
  # idx 4: 15 -> 30
  # idx 7: 2  -> 25
  res12 <- get_state_indices(1, 2, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res12, c(4L, 7L))

  # extreme (2) -> dry (0)
  # idx 5: 30 -> 0
  # idx 9: 40 -> 0
  res20 <- get_state_indices(2, 0, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res20, c(5L, 9L))
})

testthat::test_that("get_state_indices returns empty integer vector when no transitions exist", {

  prcp <- rep(0, 10)
  candidate.idx <- 1:9
  wet.thr <- 1
  extreme.thr <- 5

  res <- get_state_indices(1, 2, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res, integer(0))
})

testthat::test_that("get_state_indices correctly handles boundary values at thresholds", {

  prcp <- c(1, 1, 2, 2, 10, 10)
  candidate.idx <- 1:5
  wet.thr <- 1
  extreme.thr <- 10

  # dry -> dry
  res00 <- get_state_indices(0, 0, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res00, 1L)

  # dry -> wet
  res01 <- get_state_indices(0, 1, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res01, 2L)

  # wet -> wet
  # idx 3: 2 -> 2
  # idx 4: 2 -> 10
  # idx 5: 10 -> 10
  res11 <- get_state_indices(1, 1, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res11, c(3L, 4L, 5L))
})

testthat::test_that("get_state_indices respects candidate.idx subsetting", {

  prcp <- c(0, 5, 0, 5, 0, 5)
  wet.thr <- 1
  extreme.thr <- 10

  candidate.idx.full <- 1:5
  res.full <- get_state_indices(0, 1, prcp, candidate.idx.full, wet.thr, extreme.thr)

  candidate.idx.sub <- c(1, 3, 5)
  res.sub <- get_state_indices(0, 1, prcp, candidate.idx.sub, wet.thr, extreme.thr)

  testthat::expect_equal(res.full, c(1L, 3L, 5L))
  testthat::expect_equal(res.sub, c(1L, 2L, 3L))
})

testthat::test_that("get_state_indices output indices are relative to candidate.idx", {

  prcp <- c(0, 10, 0, 10, 0)
  candidate.idx <- c(2, 4)
  wet.thr <- 1
  extreme.thr <- 8

  res <- get_state_indices(2, 0, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res, c(1L, 2L))
})

testthat::test_that("get_state_indices works with single candidate index", {

  prcp <- c(0, 5)
  candidate.idx <- 1
  wet.thr <- 1
  extreme.thr <- 10

  res <- get_state_indices(0, 1, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_equal(res, 1L)
})

testthat::test_that("get_state_indices returns integer vector", {

  prcp <- runif(20, 0, 30)
  candidate.idx <- 1:19
  wet.thr <- 5
  extreme.thr <- 20

  res <- get_state_indices(0, 2, prcp, candidate.idx, wet.thr, extreme.thr)
  testthat::expect_type(res, "integer")
})
