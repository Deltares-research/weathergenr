
test_that("markov_next_state respects probability ordering", {

  p00 <- rep(0.6, 10)
  p01 <- rep(0.3, 10)

  p10 <- rep(0.2, 10)
  p11 <- rep(0.5, 10)

  p20 <- rep(0.1, 10)
  p21 <- rep(0.4, 10)

  # State 0
  expect_equal(markov_next_state(0, rn = 0.1, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 0)
  expect_equal(markov_next_state(0, rn = 0.7, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 1)
  expect_equal(markov_next_state(0, rn = 0.95, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 2)

  # State 1
  expect_equal(markov_next_state(1, rn = 0.1, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 0)
  expect_equal(markov_next_state(1, rn = 0.4, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 1)
  expect_equal(markov_next_state(1, rn = 0.95, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 2)

  # State 2
  expect_equal(markov_next_state(2, rn = 0.05, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 0)
  expect_equal(markov_next_state(2, rn = 0.3, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 1)
  expect_equal(markov_next_state(2, rn = 0.9, idx = 5,
                                 p00, p01, p10, p11, p20, p21), 2)
})

test_that("markov_next_state handles NA probabilities safely", {

  p00 <- rep(NA_real_, 5)
  p01 <- rep(NA_real_, 5)

  p10 <- rep(NA_real_, 5)
  p11 <- rep(NA_real_, 5)

  p20 <- rep(NA_real_, 5)
  p21 <- rep(NA_real_, 5)

  expect_silent(
    state <- markov_next_state(0, rn = 0.5, idx = 3,
                               p00, p01, p10, p11, p20, p21)
  )

  expect_true(state %in% 0:2)
})
