library(testthat)

test_that("Water year defaults to calendar year when starting month is January", {
  d <- as.Date(c("2022-01-01", "2022-12-31"))
  expect_equal(get_water_year(d, 1), c(2022, 2022))
})

test_that("Water year shifts at October (US hydrology example)", {
  d <- as.Date(c("2022-09-30", "2022-10-01", "2023-05-01"))
  # Water year starts in October: dates >= 2022-10-01 go to 2023
  expect_equal(get_water_year(d, 10), c(2022, 2023, 2023))
})

test_that("Water year shifts at June (e.g., monsoon systems)", {
  d <- as.Date(c("2022-05-31", "2022-06-01", "2023-03-01"))
  # Water year starts in June: dates >= 2022-06-01 go to 2023
  expect_equal(get_water_year(d, 6), c(2022, 2023, 2023))
})

test_that("Returns correct length and type", {
  d <- as.Date("2020-01-01") + 0:4
  wy <- get_water_year(d, 10)
  expect_length(wy, length(d))
  expect_type(wy, "integer")
})

test_that("Handles edge case for December 31st when water year is not January", {
  d <- as.Date("2022-12-31")
  expect_equal(get_water_year(d, 10), 2023)
  expect_equal(get_water_year(d, 6), 2023)
  expect_equal(get_water_year(d, 1), 2022)
})
