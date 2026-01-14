test_that("calculate_monthly_pet (hargreaves) returns correct length and no NAs for standard input", {
  month <- c(1, 6, 12)
  temp <- c(3, 18, 5)
  temp_range <- c(6, 10, 5)
  lat_deg <- 52.37

  pet <- calculate_monthly_pet(
    month = month,
    temp = temp,
    temp_range = temp_range,
    lat_deg = lat_deg,
    method = "hargreaves"
  )

  expect_length(pet, length(month))
  expect_false(anyNA(pet))
  expect_true(all(pet > 0)) # PET should be positive for these realistic inputs
})

test_that("calculate_monthly_pet (hargreaves) returns numeric vector for year-long input", {
  month <- 1:12
  temp <- seq(-2, 22, length.out = 12)
  temp_range <- rep(8, 12)
  lat_deg <- 40

  pet <- calculate_monthly_pet(
    month = month,
    temp = temp,
    temp_range = temp_range,
    lat_deg = lat_deg,
    method = "hargreaves"
  )

  expect_type(pet, "double")
  expect_length(pet, 12)
})

test_that("calculate_monthly_pet throws error for invalid months", {
  month <- c(0, 13)
  temp <- c(10, 10)
  temp_range <- c(5, 5)
  lat_deg <- 40

  expect_error(
    calculate_monthly_pet(
      month = month,
      temp = temp,
      temp_range = temp_range,
      lat_deg = lat_deg,
      method = "hargreaves"
    ),
    "'month' must contain integers between 1 and 12"
  )
})

test_that("calculate_monthly_pet throws error for invalid latitude", {
  month <- 1:3
  temp <- c(10, 15, 20)
  temp_range <- c(6, 7, 8)

  expect_error(
    calculate_monthly_pet(month, temp, temp_range, -91, method = "hargreaves"),
    "'lat_deg' must be a finite numeric value between -90 and 90"
  )
  expect_error(
    calculate_monthly_pet(month, temp, temp_range, 100, method = "hargreaves"),
    "'lat_deg' must be a finite numeric value between -90 and 90"
  )
  expect_error(
    calculate_monthly_pet(month, temp, temp_range, "north", method = "hargreaves"),
    "'lat_deg' must be a finite numeric value between -90 and 90"
  )
})

test_that("calculate_monthly_pet (hargreaves) returns values for Southern Hemisphere", {
  month <- c(1, 7, 12)
  temp <- c(25, 10, 27)
  temp_range <- c(7, 4, 6)
  lat_deg <- -33.9 # e.g., Sydney

  pet <- calculate_monthly_pet(
    month = month,
    temp = temp,
    temp_range = temp_range,
    lat_deg = lat_deg,
    method = "hargreaves"
  )

  expect_length(pet, 3)
  expect_true(all(pet > 0))
})

test_that("calculate_monthly_pet coerces month to integers", {
  out1 <- calculate_monthly_pet(
    month = c(1.5, 2),
    temp = c(10, 10),
    temp_range = c(5, 5),
    lat_deg = 40,
    method = "hargreaves"
  )

  out2 <- calculate_monthly_pet(
    month = c(1, 2),
    temp = c(10, 10),
    temp_range = c(5, 5),
    lat_deg = 40,
    method = "hargreaves"
  )

  expect_equal(out1, out2)
})
