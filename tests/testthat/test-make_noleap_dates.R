test_that("returns correct length and Date class", {

  out <- make_noleap_dates("1980-01-01", 365)

  expect_s3_class(out, "Date")
  expect_length(out, 365)
})


test_that("never includes February 29", {

  out <- make_noleap_dates("1980-01-01", 2000)

  expect_false(any(format(out, "%m-%d") == "02-29"))
})


test_that("start_date on Feb 29 is normalized to Feb 28", {

  out <- make_noleap_dates("2000-02-29", 3)

  expect_identical(
    format(out[1], "%m-%d"),
    "02-28"
  )
})


test_that("dates progress monotonically without gaps", {

  out <- make_noleap_dates("1999-12-30", 10)

  # Differences should always be 1 day
  expect_true(all(diff(out) == 1))
})


test_that("year increments exactly every 365 days", {

  start <- as.Date("2010-03-01")
  out   <- make_noleap_dates(start, 800)

  years <- as.integer(format(out, "%Y"))

  # Year should increase only when synthetic DOY wraps
  year_jumps <- which(diff(years) > 0)

  expect_true(all(diff(year_jumps) >= 364))
  expect_true(all(diff(year_jumps) <= 366))
})


test_that("month-day sequence repeats identically every 365 days", {

  out <- make_noleap_dates("2015-01-01", 1000)

  md <- format(out, "%m-%d")

  expect_identical(
    md[1:365],
    md[366:730]
  )
})


test_that("start date aligns correctly with template", {

  out <- make_noleap_dates("2021-06-15", 5)

  expect_identical(
    format(out, "%m-%d"),
    c("06-15", "06-16", "06-17", "06-18", "06-19")
  )
})


test_that("supports n_dates = 1 edge case", {

  out <- make_noleap_dates("1990-01-01", 1)

  expect_length(out, 1)
  expect_identical(out, as.Date("1990-01-01"))
})


test_that("long sequences remain consistent over multiple synthetic years", {

  out <- make_noleap_dates("1985-07-01", 5000)

  # Still no leap days in month-day labels
  expect_false(any(format(out, "%m-%d") == "02-29"))

  # Dates must be strictly increasing
  expect_true(all(diff(out) > 0))

  # Differences should be 1 day, except 2 days when skipping Feb 29 in leap years
  d <- as.integer(diff(out))
  expect_true(all(d %in% c(1L, 2L)))

  # All "2-day jumps" must be exactly Feb 28 -> Mar 01 in a leap year
  jump2_idx <- which(d == 2L)
  if (length(jump2_idx) > 0) {
    prev_dates <- out[jump2_idx]
    next_dates <- out[jump2_idx + 1L]

    expect_true(all(format(prev_dates, "%m-%d") == "02-28"))
    expect_true(all(format(next_dates, "%m-%d") == "03-01"))

    prev_years <- as.integer(format(prev_dates, "%Y"))
    is_leap <- (prev_years %% 4 == 0 & prev_years %% 100 != 0) | (prev_years %% 400 == 0)
    expect_true(all(is_leap))
  }
})

