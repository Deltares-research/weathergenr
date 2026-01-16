# Tests for match_transition_positions function
# This function returns positions within day0_idx where a specified state transition occurs

testthat::test_that("match_transition_positions returns correct transitions for simple example", {
  prcp <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
  day0_idx <- 1:(length(prcp) - 1)
  wet_threshold <- 1
  extreme_threshold <- 20
  # dry (0) -> wet (1)
  # transitions at:
  # idx 2: 0 -> 5
  # idx 6: 0 -> 2
  res01 <- match_transition_positions(0, 1, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res01, c(2L, 6L))
  # wet (1) -> extreme (2)
  # idx 4: 15 -> 30
  # idx 7: 2  -> 25
  res12 <- match_transition_positions(1, 2, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res12, c(4L, 7L))
  # extreme (2) -> dry (0)
  # idx 5: 30 -> 0
  # idx 9: 40 -> 0
  res20 <- match_transition_positions(2, 0, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res20, c(5L, 9L))
})

testthat::test_that("match_transition_positions returns empty integer vector when no transitions exist", {
  prcp <- rep(0, 10)
  day0_idx <- 1:9
  wet_threshold <- 1
  extreme_threshold <- 5
  res <- match_transition_positions(1, 2, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res, integer(0))
})

testthat::test_that("match_transition_positions correctly handles boundary values at thresholds", {
  prcp <- c(1, 1, 2, 2, 10, 10)
  day0_idx <- 1:5
  wet_threshold <- 1
  extreme_threshold <- 10
  # dry -> dry
  res00 <- match_transition_positions(0, 0, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res00, 1L)
  # dry -> wet
  res01 <- match_transition_positions(0, 1, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res01, 2L)
  # wet -> wet
  # idx 3: 2 -> 2
  # idx 4: 2 -> 10
  # idx 5: 10 -> 10
  res11 <- match_transition_positions(1, 1, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res11, c(3L, 4L, 5L))
})

testthat::test_that("match_transition_positions respects day0_idx subsetting", {
  prcp <- c(0, 5, 0, 5, 0, 5)
  wet_threshold <- 1
  extreme_threshold <- 10
  day0_idx_full <- 1:5
  res_full <- match_transition_positions(0, 1, prcp, day0_idx_full, wet_threshold, extreme_threshold)
  day0_idx_sub <- c(1, 3, 5)
  res_sub <- match_transition_positions(0, 1, prcp, day0_idx_sub, wet_threshold, extreme_threshold)
  testthat::expect_equal(res_full, c(1L, 3L, 5L))
  testthat::expect_equal(res_sub, c(1L, 2L, 3L))
})

testthat::test_that("match_transition_positions output indices are relative to day0_idx", {
  prcp <- c(0, 10, 0, 10, 0)
  day0_idx <- c(2, 4)
  wet_threshold <- 1
  extreme_threshold <- 8
  res <- match_transition_positions(2, 0, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res, c(1L, 2L))
})

testthat::test_that("match_transition_positions works with single candidate index", {
  prcp <- c(0, 5)
  day0_idx <- 1
  wet_threshold <- 1
  extreme_threshold <- 10
  res <- match_transition_positions(0, 1, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_equal(res, 1L)
})

testthat::test_that("match_transition_positions returns integer vector", {
  prcp <- runif(20, 0, 30)
  day0_idx <- 1:19
  wet_threshold <- 5
  extreme_threshold <- 20
  res <- match_transition_positions(0, 2, prcp, day0_idx, wet_threshold, extreme_threshold)
  testthat::expect_type(res, "integer")
})
