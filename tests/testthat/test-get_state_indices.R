test_that("get_state_indices returns correct transitions", {

  PRCP <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
  cur_day <- 1:(length(PRCP) - 1)

  wet_thr <- 1
  extreme_thr <- 20

  # Dry -> Dry
  idx <- get_state_indices(0, 0, PRCP, cur_day, wet_thr, extreme_thr)
  expect_true(all(PRCP[cur_day[idx]] <= wet_thr))
  expect_true(all(PRCP[cur_day[idx] + 1] <= wet_thr))

  # Dry -> Wet
  idx <- get_state_indices(0, 1, PRCP, cur_day, wet_thr, extreme_thr)
  expect_true(all(PRCP[cur_day[idx]] <= wet_thr))
  expect_true(all(PRCP[cur_day[idx] + 1] > wet_thr &
                    PRCP[cur_day[idx] + 1] <= extreme_thr))

  # Wet -> Extreme
  idx <- get_state_indices(1, 2, PRCP, cur_day, wet_thr, extreme_thr)
  expect_true(all(PRCP[cur_day[idx]] > wet_thr &
                    PRCP[cur_day[idx]] <= extreme_thr))
  expect_true(all(PRCP[cur_day[idx] + 1] > extreme_thr))

  # Extreme -> Dry
  idx <- get_state_indices(2, 0, PRCP, cur_day, wet_thr, extreme_thr)
  expect_true(all(PRCP[cur_day[idx]] > extreme_thr))
  expect_true(all(PRCP[cur_day[idx] + 1] <= wet_thr))
})

test_that("get_state_indices handles empty results safely", {

  PRCP <- c(0, 0, 0, 0)
  cur_day <- 1:(length(PRCP) - 1)

  idx <- get_state_indices(2, 2, PRCP, cur_day, wet_thr = 1, extreme_thr = 10)
  expect_length(idx, 0)
})

