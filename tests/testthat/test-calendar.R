# ==============================================================================
# Unit Tests for calendar.R
# ==============================================================================
# Tests for: generate_noleap_dates, compute_water_year, find_leap_day_indices,
#            normalize_calendar, build_historical_dates
# ==============================================================================

library(testthat)

# ==============================================================================
# Tests for generate_noleap_dates()
# ==============================================================================

test_that("generate_noleap_dates returns correct length", {
  result <- generate_noleap_dates("1980-01-01", 365)
  expect_length(result, 365)

  result <- generate_noleap_dates("1980-01-01", 1000)
  expect_length(result, 1000)

  result <- generate_noleap_dates("2000-06-15", 1)
  expect_length(result, 1)
})

test_that("generate_noleap_dates returns Date class", {
            result <- generate_noleap_dates("1980-01-01", 100)
            expect_s3_class(result, "Date")
})

test_that("generate_noleap_dates excludes February 29", {

  # Generate 3 full years (1095 days)
  result <- generate_noleap_dates("1980-01-01", 1095)
  month_days <- format(result, "%m-%d")
  expect_false("02-29" %in% month_days)
})

test_that("generate_noleap_dates starts from correct date", {
  result <- generate_noleap_dates("1995-03-15", 10)
  expect_equal(result[1], as.Date("1995-03-15"))

  result <- generate_noleap_dates("2020-12-31", 5)
  expect_equal(result[1], as.Date("2020-12-31"))
})

test_that("generate_noleap_dates handles Feb 29 start date by mapping to Feb 28", {
  # Feb 29 should be mapped to Feb 28

  result <- generate_noleap_dates("2000-02-29", 5)
  expect_equal(format(result[1], "%m-%d"), "02-28")
})

test_that("generate_noleap_dates increments year correctly", {
  # Start from Dec 31, next day should be Jan 1 of next year
  result <- generate_noleap_dates("1999-12-31", 3)
  expect_equal(result[1], as.Date("1999-12-31"))
  expect_equal(result[2], as.Date("2000-01-01"))
  expect_equal(result[3], as.Date("2000-01-02"))
})

test_that("generate_noleap_dates produces 365 unique month-day combinations per year", {
  result <- generate_noleap_dates("2001-01-01", 365)
  month_days <- format(result, "%m-%d")
  expect_equal(length(unique(month_days)), 365)
})

test_that("generate_noleap_dates accepts character date input", {
  result <- generate_noleap_dates("1980-01-01", 10)
  expect_s3_class(result, "Date")
  expect_length(result, 10)
})

test_that("generate_noleap_dates accepts Date object input", {
  result <- generate_noleap_dates(as.Date("1980-01-01"), 10)
  expect_s3_class(result, "Date")
  expect_length(result, 10)
})

test_that("generate_noleap_dates errors on invalid start_date", {
  expect_error(
    generate_noleap_dates(NA, 10),
    "must be a valid Date"
  )
  expect_error(
    generate_noleap_dates("invalid-date", 10),
    "must be a valid Date"
  )
  expect_error(
    generate_noleap_dates(c("1980-01-01", "1980-01-02"), 10),
    "must be a valid Date"
  )
  expect_error(
    generate_noleap_dates(NULL, 10),
    "must be a valid Date"
  )
})

test_that("generate_noleap_dates errors on invalid n_days", {
  expect_error(
    generate_noleap_dates("1980-01-01", 0),
    "must be a positive integer"
  )
  expect_error(
    generate_noleap_dates("1980-01-01", -5),
    "must be a positive integer"
  )
  expect_error(
    generate_noleap_dates("1980-01-01", 3.5),
    "must be a positive integer"
  )
  expect_error(
    generate_noleap_dates("1980-01-01", "ten"),
    "must be a positive integer"
  )
  expect_error(
    generate_noleap_dates("1980-01-01", NA),
    "must be a positive integer"
  )
  expect_error(
    generate_noleap_dates("1980-01-01", Inf),
    "must be a positive integer"
  )
})

test_that("generate_noleap_dates produces continuous sequence", {
  result <- generate_noleap_dates("2001-01-01", 365)
  # Each day should be exactly 1 day after the previous (in no-leap calendar)
  # Check that we don't skip any month-day combinations
  expected_md <- format(
    seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day"),
    "%m-%d"
  )
  expected_md <- expected_md[expected_md != "02-29"]
  result_md <- format(result, "%m-%d")
  expect_equal(result_md, expected_md)
})


# ==============================================================================
# Tests for compute_water_year()
# ==============================================================================

test_that("compute_water_year returns integer vector", {
  dates <- as.Date(c("2022-01-01", "2022-06-15", "2022-12-31"))
  result <- compute_water_year(dates, water_year_start_month = 1)
  expect_type(result, "integer")
  expect_length(result, 3)
})

test_that("compute_water_year with January start equals calendar year", {
  dates <- as.Date(c("2020-01-01", "2020-06-15", "2020-12-31", "2021-01-01"))
  result <- compute_water_year(dates, water_year_start_month = 1)
  expected <- c(2020L, 2020L, 2020L, 2021L)
  expect_equal(result, expected)
})

test_that("compute_water_year with October start works correctly", {
  # Standard water year: Oct 1 - Sep 30
  dates <- as.Date(c("2022-09-30", "2022-10-01", "2023-01-15", "2023-09-30"))
  result <- compute_water_year(dates, water_year_start_month = 10)
  # Sep 30 2022 belongs to WY 2022
  # Oct 1 2022 belongs to WY 2023
  # Jan 15 2023 belongs to WY 2023
  # Sep 30 2023 belongs to WY 2023

  expected <- c(2022L, 2023L, 2023L, 2023L)
  expect_equal(result, expected)
})

test_that("compute_water_year with June start works correctly", {
  dates <- as.Date(c("2022-05-31", "2022-06-01", "2022-12-15", "2023-05-31"))
  result <- compute_water_year(dates, water_year_start_month = 6)
  # May 31 2022 belongs to WY 2022
  # Jun 1 2022 belongs to WY 2023
  # Dec 15 2022 belongs to WY 2023
  # May 31 2023 belongs to WY 2023
  expected <- c(2022L, 2023L, 2023L, 2023L)
  expect_equal(result, expected)
})

test_that("compute_water_year handles boundary dates", {
  # Test with December start (month 12)
  dates <- as.Date(c("2022-11-30", "2022-12-01", "2023-11-30"))
  result <- compute_water_year(dates, water_year_start_month = 12)
  expected <- c(2022L, 2023L, 2023L)
  expect_equal(result, expected)
})

test_that("compute_water_year handles single date", {
  result <- compute_water_year(as.Date("2022-10-15"), water_year_start_month = 10)
  expect_equal(result, 2023L)
})

test_that("compute_water_year handles empty date vector", {
  result <- compute_water_year(as.Date(character(0)), water_year_start_month = 1)
  expect_length(result, 0)
  expect_type(result, "integer")
})

test_that("compute_water_year errors on non-Date input", {
  expect_error(
    compute_water_year("2022-01-01", water_year_start_month = 1),
    "must be of class Date"
  )
  expect_error(
    compute_water_year(c(2022, 2023), water_year_start_month = 1),
    "must be of class Date"
  )
})

test_that("compute_water_year errors on invalid water_year_start_month", {
  dates <- as.Date("2022-01-01")

  expect_error(
    compute_water_year(dates, water_year_start_month = 0),
    "must be an integer between 1 and 12"
  )
  expect_error(
    compute_water_year(dates, water_year_start_month = 13),
    "must be an integer between 1 and 12"
  )
  expect_error(
    compute_water_year(dates, water_year_start_month = -1),
    "must be an integer between 1 and 12"
  )
  expect_error(
    compute_water_year(dates, water_year_start_month = 6.5),
    "must be an integer between 1 and 12"
  )
  expect_error(
    compute_water_year(dates, water_year_start_month = NA),
    "must be an integer between 1 and 12"
  )
  expect_error(
    compute_water_year(dates, water_year_start_month = "October"),
    "must be an integer between 1 and 12"
  )
})

test_that("compute_water_year default is calendar year", {
  dates <- as.Date(c("2022-06-15", "2022-12-31"))
  result <- compute_water_year(dates)
  expect_equal(result, c(2022L, 2022L))
})


# ==============================================================================
# Tests for find_leap_day_indices()
# ==============================================================================

test_that("find_leap_day_indices finds leap days", {
  dates <- as.Date(c("1980-02-28", "1980-02-29", "1980-03-01"))
  result <- find_leap_day_indices(dates)
  expect_equal(result, 2L)
})

test_that("find_leap_day_indices finds multiple leap days", {
  dates <- as.Date(c(
    "2000-02-29", "2000-03-01",
    "2004-02-29", "2004-03-01",
    "2008-02-29"
  ))
  result <- find_leap_day_indices(dates)
  expect_equal(result, c(1L, 3L, 5L))
})

test_that("find_leap_day_indices returns NULL when no leap days", {
  dates <- as.Date(c("2001-01-01", "2001-02-28", "2001-12-31"))
  result <- find_leap_day_indices(dates)
  expect_null(result)
})

test_that("find_leap_day_indices handles single leap day date", {
  result <- find_leap_day_indices(as.Date("2000-02-29"))
  expect_equal(result, 1L)
})

test_that("find_leap_day_indices handles single non-leap day date", {
  result <- find_leap_day_indices(as.Date("2000-02-28"))
  expect_null(result)
})

test_that("find_leap_day_indices returns integer vector", {
  dates <- as.Date(c("2000-02-29", "2004-02-29"))
  result <- find_leap_day_indices(dates)
  expect_type(result, "integer")
})

test_that("find_leap_day_indices accepts character input", {
  result <- find_leap_day_indices(c("2000-02-28", "2000-02-29", "2000-03-01"))
  expect_equal(result, 2L)
})

test_that("find_leap_day_indices errors on NA values", {
  expect_error(
    find_leap_day_indices(c(as.Date("2000-02-29"), NA)),
    "must be coercible to Date with no NA values"
  )
})

test_that("find_leap_day_indices errors on non-coercible input", {
  expect_error(
    find_leap_day_indices(c("not-a-date", "also-not")),
    "must be coercible to Date with no NA values"
  )
})


# ==============================================================================
# Tests for normalize_calendar()
# ==============================================================================

test_that("normalize_calendar removes leap days from dates", {
  dates <- as.Date(c("2000-02-28", "2000-02-29", "2000-03-01"))
  data <- list(
    cell1 = data.frame(precip = c(1.0, 2.0, 3.0), temp = c(10, 11, 12))
  )

  result <- normalize_calendar(dates, data, verbose = FALSE)

  expect_length(result$dates, 2)
  expect_equal(result$dates, as.Date(c("2000-02-28", "2000-03-01")))
  expect_equal(nrow(result$data$cell1), 2)
  expect_equal(result$data$cell1$precip, c(1.0, 3.0))
})

test_that("normalize_calendar handles multiple grid cells", {
  dates <- as.Date(c("2000-02-28", "2000-02-29", "2000-03-01"))
  data <- list(
    cell1 = data.frame(precip = c(1, 2, 3)),
    cell2 = data.frame(precip = c(4, 5, 6)),
    cell3 = data.frame(precip = c(7, 8, 9))
  )

  result <- normalize_calendar(dates, data, verbose = FALSE)

  expect_equal(names(result$data), c("cell1", "cell2", "cell3"))
  expect_equal(nrow(result$data$cell1), 2)
  expect_equal(nrow(result$data$cell2), 2)
  expect_equal(nrow(result$data$cell3), 2)
  expect_equal(result$data$cell1$precip, c(1, 3))
  expect_equal(result$data$cell2$precip, c(4, 6))
})

test_that("normalize_calendar returns unchanged data when no leap days", {
  dates <- as.Date(c("2001-02-27", "2001-02-28", "2001-03-01"))
  data <- list(
    cell1 = data.frame(precip = c(1, 2, 3))
  )

  result <- normalize_calendar(dates, data, verbose = FALSE)

  expect_equal(result$dates, dates)
  expect_equal(result$data$cell1$precip, c(1, 2, 3))
})

test_that("normalize_calendar removes multiple leap days", {
  dates <- as.Date(c(
    "2000-02-28", "2000-02-29", "2000-03-01",
    "2004-02-28", "2004-02-29", "2004-03-01"
  ))
  data <- list(
    cell1 = data.frame(val = 1:6)
  )

  result <- normalize_calendar(dates, data, verbose = FALSE)

  expect_length(result$dates, 4)
  expect_false("02-29" %in% format(result$dates, "%m-%d"))
  expect_equal(result$data$cell1$val, c(1, 3, 4, 6))
})

test_that("normalize_calendar preserves column structure", {
  dates <- as.Date(c("2000-02-28", "2000-02-29", "2000-03-01"))
  data <- list(
    cell1 = data.frame(
      precip = c(1, 2, 3),
      tmax = c(10, 11, 12),
      tmin = c(0, 1, 2),
      wind = c(5, 6, 7)
    )
  )

  result <- normalize_calendar(dates, data, verbose = FALSE)

  expect_equal(names(result$data$cell1), c("precip", "tmax", "tmin", "wind"))
})

test_that("normalize_calendar errors on non-data.frame elements", {
  dates <- as.Date(c("2000-02-28", "2000-02-29", "2000-03-01"))
  data <- list(
    cell1 = c(1, 2, 3)  # vector, not data.frame
  )

  expect_error(
    normalize_calendar(dates, data, verbose = FALSE),
    "must be a data.frame"
  )
})

test_that("normalize_calendar errors on mismatched row counts", {
  dates <- as.Date(c("2000-02-28", "2000-02-29", "2000-03-01"))
  data <- list(
    cell1 = data.frame(precip = c(1, 2))  # only 2 rows, not 3
  )

  expect_error(
    normalize_calendar(dates, data, verbose = FALSE),
    "must have nrow equal to"
  )
})

test_that("normalize_calendar returns correct list structure", {
  dates <- as.Date(c("2000-02-28", "2000-03-01"))
  data <- list(cell1 = data.frame(val = c(1, 2)))

  result <- normalize_calendar(dates, data, verbose = FALSE)

  expect_type(result, "list")
  expect_named(result, c("dates", "data"))
  expect_s3_class(result$dates, "Date")
  expect_type(result$data, "list")
})


# ==============================================================================
# Tests for build_historical_dates()
# ==============================================================================

test_that("build_historical_dates returns correct structure", {
  # Create 2 complete years of data (730 days, no leap days)
  dates <- generate_noleap_dates("2001-01-01", 730)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  expect_type(result, "list")
  expect_named(result, c("dates_df", "wyear_idx", "complete_wyears"))
  expect_s3_class(result$dates_df, "tbl_df")
  expect_type(result$wyear_idx, "integer")
  expect_type(result$complete_wyears, "integer")
})

test_that("build_historical_dates identifies complete years", {
  # Create exactly 3 complete years
  dates <- generate_noleap_dates("2001-01-01", 365 * 3)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  expect_length(result$complete_wyears, 3)
  expect_equal(result$complete_wyears, c(2001L, 2002L, 2003L))
})

test_that("build_historical_dates excludes incomplete years", {
  # Create 2.5 years (incomplete third year)
  dates <- generate_noleap_dates("2001-01-01", 365 * 2 + 180)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  # Should only have 2 complete years
  expect_length(result$complete_wyears, 2)
  expect_equal(nrow(result$dates_df), 730)
})

test_that("build_historical_dates works with October water year", {
  # Create dates spanning Oct 2000 to Sep 2003 (3 complete water years)
  # Start from Oct 1, 2000
  dates <- generate_noleap_dates("2000-10-01", 365 * 3)

  result <- build_historical_dates(dates, year_start_month = 10, verbose = FALSE)

  expect_length(result$complete_wyears, 3)
  # Water years should be 2001, 2002, 2003 (Oct-Sep pattern)
  expect_equal(result$complete_wyears, c(2001L, 2002L, 2003L))
})

test_that("build_historical_dates dates_df has correct columns", {
  dates <- generate_noleap_dates("2001-01-01", 365)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  expect_named(
    result$dates_df,
    c("dateo", "year", "wyear", "month", "day", "date")
  )
})

test_that("build_historical_dates wyear_idx maps correctly", {
  dates <- generate_noleap_dates("2001-01-01", 365)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  # wyear_idx should map back to original dates
  expect_equal(dates[result$wyear_idx], result$dates_df$dateo)
})

test_that("build_historical_dates finds longest contiguous block", {
  # Create dates with a gap (simulating incomplete year in middle)
  # 2001: complete, 2002: incomplete (missing), 2003-2005: complete
  dates_2001 <- generate_noleap_dates("2001-01-01", 365)
  dates_2002_partial <- generate_noleap_dates("2002-01-01", 100)  # incomplete
  dates_2003_2005 <- generate_noleap_dates("2003-01-01", 365 * 3)

  dates <- c(dates_2001, dates_2002_partial, dates_2003_2005)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  # Should select longest contiguous block: 2003-2005 (3 years)
  expect_length(result$complete_wyears, 3)
  expect_equal(result$complete_wyears, c(2003L, 2004L, 2005L))
})

test_that("build_historical_dates errors with no complete years", {
  # Only 100 days - less than one complete year
  dates <- generate_noleap_dates("2001-01-01", 100)

  expect_error(
    build_historical_dates(dates, year_start_month = 1, verbose = FALSE),
    "No complete years found"
  )
})

test_that("build_historical_dates handles June water year start", {
  # Start from June 1, create 2 complete water years
  dates <- generate_noleap_dates("2001-06-01", 365 * 2)

  result <- build_historical_dates(dates, year_start_month = 6, verbose = FALSE)

  expect_length(result$complete_wyears, 2)
  # June 2001 starts WY 2002
  expect_equal(result$complete_wyears, c(2002L, 2003L))
})

test_that("build_historical_dates dateo matches original dates", {
  dates <- generate_noleap_dates("2001-01-01", 365)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  # All dateo values should be present in original dates

  expect_true(all(result$dates_df$dateo %in% dates))
})

test_that("build_historical_dates month and day columns are correct", {
  dates <- generate_noleap_dates("2001-01-01", 365)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  # Check first and last entries
  expect_equal(result$dates_df$month[1], 1L)
  expect_equal(result$dates_df$day[1], 1L)
  expect_equal(result$dates_df$month[365], 12L)
  expect_equal(result$dates_df$day[365], 31L)
})

test_that("build_historical_dates year column is calendar year", {
  dates <- generate_noleap_dates("2001-01-01", 365)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  expect_true(all(result$dates_df$year == 2001L))
})


# ==============================================================================
# Integration Tests
# ==============================================================================

test_that("normalize_calendar and build_historical_dates work together", {
  # Create dates including leap days
  dates <- seq.Date(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "day")
  data <- list(
    cell1 = data.frame(val = seq_along(dates))
  )

  # First normalize (remove leap days)
  normalized <- normalize_calendar(dates, data, verbose = FALSE)

  # Then build historical dates
  result <- build_historical_dates(
    normalized$dates,
    year_start_month = 1,
    verbose = FALSE
  )

  # Should have 3 complete years
  expect_length(result$complete_wyears, 3)
  expect_equal(result$complete_wyears, c(2000L, 2001L, 2002L))
})

test_that("generate_noleap_dates output works with build_historical_dates", {
  dates <- generate_noleap_dates("2001-01-01", 365 * 5)

  result <- build_historical_dates(dates, year_start_month = 1, verbose = FALSE)

  expect_length(result$complete_wyears, 5)
})

test_that("compute_water_year consistent with build_historical_dates", {
  dates <- generate_noleap_dates("2000-10-01", 365 * 2)

  # Direct water year computation
  direct_wyear <- compute_water_year(dates, water_year_start_month = 10)

  # Via build_historical_dates
  result <- build_historical_dates(dates, year_start_month = 10, verbose = FALSE)

  # The wyear values in dates_df should match direct computation
  idx <- result$wyear_idx
  expect_equal(result$dates_df$wyear, direct_wyear[idx])
})


# ==============================================================================
# Edge Case Tests
# ==============================================================================

test_that("functions handle year 2000 (leap year) correctly", {
  # 2000 was a leap year (divisible by 400)
  result <- find_leap_day_indices(as.Date("2000-02-29"))
  expect_equal(result, 1L)
})

test_that("functions handle year 1900 (not leap year) correctly", {
  # 1900 was NOT a leap year (divisible by 100 but not 400)
  dates <- seq.Date(as.Date("1900-02-27"), as.Date("1900-03-02"), by = "day")
  result <- find_leap_day_indices(dates)
  expect_null(result)  # No Feb 29 in 1900
})

test_that("functions handle year 2100 (not leap year) correctly", {
  # 2100 will NOT be a leap year
  dates <- seq.Date(as.Date("2100-02-27"), as.Date("2100-03-02"), by = "day")
  result <- find_leap_day_indices(dates)
  expect_null(result)
})

test_that("compute_water_year handles dates across multiple decades", {
  dates <- as.Date(c("1990-10-15", "2000-10-15", "2010-10-15", "2020-10-15"))
  result <- compute_water_year(dates, water_year_start_month = 10)
  expect_equal(result, c(1991L, 2001L, 2011L, 2021L))
})

test_that("generate_noleap_dates handles very long sequences", {
  # 100 years of data
  result <- generate_noleap_dates("1900-01-01", 365 * 100)
  expect_length(result, 36500)

  # Should never contain Feb 29
  expect_false("02-29" %in% format(result, "%m-%d"))
})
