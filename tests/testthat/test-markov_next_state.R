test_that("markov_next_state respects probability ordering", {

  p00 <- rep(0.6, 10)
  p01 <- rep(0.3, 10)

  p10 <- rep(0.2, 10)
  p11 <- rep(0.5, 10)

  p20 <- rep(0.1, 10)
  p21 <- rep(0.4, 10)

  # State 0
  expect_equal(weathergenr::markov_next_state(0, u_rand = 0.1, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 0)
  expect_equal(weathergenr::markov_next_state(0, u_rand = 0.7, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 1)
  expect_equal(weathergenr::markov_next_state(0, u_rand = 0.95, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 2)

  # State 1
  expect_equal(weathergenr::markov_next_state(1, u_rand = 0.1, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 0)
  expect_equal(weathergenr::markov_next_state(1, u_rand = 0.4, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 1)
  expect_equal(weathergenr::markov_next_state(1, u_rand = 0.95, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 2)

  # State 2
  expect_equal(weathergenr::markov_next_state(2, u_rand = 0.05, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 0)
  expect_equal(weathergenr::markov_next_state(2, u_rand = 0.3, idx = 5,
                                              p00, p01, p10, p11, p20, p21), 1)
  expect_equal(weathergenr::markov_next_state(2, u_rand = 0.9, idx = 5,
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
    state <- weathergenr::markov_next_state(0, u_rand = 0.5, idx = 3,
                                            p00, p01, p10, p11, p20, p21)
  )

  expect_true(state %in% 0:2)
})

test_that("markov_next_state clamps probabilities when sum exceeds 1", {

  # Invalid probabilities (sum > 1)
  p00 <- rep(0.8, 5)
  p01 <- rep(0.6, 5)  # sum = 1.4
  p10 <- rep(0.8, 5)
  p11 <- rep(0.6, 5)
  p20 <- rep(0.8, 5)
  p21 <- rep(0.6, 5)

  # Suppress the repeated clamp wau_randing spam from replicate()
  out <- suppressWarnings(
    replicate(100, weathergenr::markov_next_state(
      state_prev = 0,
      u_rand = runif(1),
      idx = 3,
      p00, p01, p10, p11, p20, p21
    ))
  )

  expect_true(all(out %in% 0:2))
})

test_that("markov_next_state handles invalid idx gracefully", {

  p00 <- rep(0.6, 3)
  p01 <- rep(0.2, 3)
  p10 <- rep(0.3, 3)
  p11 <- rep(0.4, 3)
  p20 <- rep(0.1, 3)
  p21 <- rep(0.3, 3)

  # idx too small
  out1 <- suppressWarnings(
    weathergenr::markov_next_state(0, 0.5, idx = -5, p00, p01, p10, p11, p20, p21)
  )
  expect_true(out1 %in% 0:2)

  # idx too large
  out2 <- suppressWarnings(
    weathergenr::markov_next_state(1, 0.5, idx = 99, p00, p01, p10, p11, p20, p21)
  )
  expect_true(out2 %in% 0:2)
})

test_that("markov_next_state handles invalid state_prev safely", {

  p00 <- rep(0.6, 5)
  p01 <- rep(0.2, 5)
  p10 <- rep(0.3, 5)
  p11 <- rep(0.4, 5)
  p20 <- rep(0.1, 5)
  p21 <- rep(0.3, 5)

  out <- suppressWarnings(
    weathergenr::markov_next_state(
      state_prev = 99,  # invalid
      u_rand = 0.5,
      idx = 2,
      p00, p01, p10, p11, p20, p21
    )
  )

  expect_true(out %in% 0:2)
})

test_that("markov_next_state is reproducible for fixed u_rand", {

  p00 <- rep(0.6, 10)
  p01 <- rep(0.2, 10)
  p10 <- rep(0.3, 10)
  p11 <- rep(0.4, 10)
  p20 <- rep(0.1, 10)
  p21 <- rep(0.3, 10)

  u_rand <- 0.37

  out1 <- weathergenr::markov_next_state(0, u_rand, 5, p00, p01, p10, p11, p20, p21)
  out2 <- weathergenr::markov_next_state(0, u_rand, 5, p00, p01, p10, p11, p20, p21)

  expect_identical(out1, out2)
})
