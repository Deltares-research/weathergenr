library(testthat)

test_that("pet_hargreaves returns correct length and no NAs for standard input", {
  months <- c(1, 6, 12)
  temp <- c(3, 18, 5)
  tdiff <- c(6, 10, 5)
  lat <- 52.37
  pet <- pet_hargreaves(months, temp, tdiff, lat)
  expect_length(pet, length(months))
  expect_false(anyNA(pet))
  expect_true(all(pet > 0)) # PET should be positive for these realistic inputs
})

test_that("pet_hargreaves returns numeric vector for year-long input", {
  months <- 1:12
  temp <- seq(-2, 22, length.out = 12)
  tdiff <- rep(8, 12)
  lat <- 40
  pet <- pet_hargreaves(months, temp, tdiff, lat)
  expect_type(pet, "double")
  expect_length(pet, 12)
})

test_that("pet_hargreaves throws error for invalid months", {
  months <- c(0, 13)
  temp <- c(10, 10)
  tdiff <- c(5, 5)
  lat <- 40
  expect_error(
    pet_hargreaves(months, temp, tdiff, lat),
    "months.*1 and 12"
  )
})

test_that("pet_hargreaves throws error for invalid latitude", {
  months <- 1:3
  temp <- c(10, 15, 20)
  tdiff <- c(6, 7, 8)
  expect_error(pet_hargreaves(months, temp, tdiff, -91), "lat.*-90 and 90")
  expect_error(pet_hargreaves(months, temp, tdiff, 100), "lat.*-90 and 90")
  expect_error(pet_hargreaves(months, temp, tdiff, "north"), "lat.*numeric")
})

test_that("pet_hargreaves returns values for Southern Hemisphere", {
  months <- c(1, 7, 12)
  temp <- c(25, 10, 27)
  tdiff <- c(7, 4, 6)
  lat <- -33.9 # e.g., Sydney
  pet <- pet_hargreaves(months, temp, tdiff, lat)
  expect_length(pet, 3)
  expect_true(all(pet > 0))
})

test_that("pet_hargreaves coerces months to integers", {
  out1 <- pet_hargreaves(c(1.5, 2), c(10, 10), c(5, 5), 40)
  out2 <- pet_hargreaves(c(1,   2), c(10, 10), c(5, 5), 40)
  expect_equal(out1, out2)
})
