#' @title Generate 365-Day Yearly Date Sequence (No Leap Days)
#'
#' @description
#' Generates a sequence of daily dates for a specified number of years, always omitting February 29 (leap days).
#' This utility is especially useful for stochastic weather generators and synthetic climate simulations where every year
#' must have exactly 365 days, regardless of actual leap years in the calendar.
#'
#' @param start_year Integer. The starting year of the date sequence (e.g., 2000).
#' @param n_years Integer. The number of years to generate (e.g., 10).
#'
#' @return
#' A vector of `Date` objects, length `365 * n_years`, with each year containing exactly 365 days (no Feb 29).
#'
#' @examples
#' # Generate 2 synthetic years starting from 2000, with no leap days
#' dates <- make_noleap_dates(2000, 2)
#' head(dates)
#' tail(dates)
#' any(format(dates, "%m-%d") == "02-29") # should be FALSE
#' length(dates) # should be 730 (2 years * 365)
#'
#' @export
make_noleap_dates <- function(start_year, n_years) {
  # Build a template year (non-leap, e.g. 2001) and remove Feb 29 if present
  template <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  template <- template[format(template, "%m-%d") != "02-29"]  # Removes Feb 29 if present (should not be, but robust)

  dates <- as.Date(unlist(lapply(0:(n_years - 1), function(i) {
    # Use template month and day, but with incremented year
    sprintf("%04d-%02d-%02d",
            start_year + i,
            as.integer(format(template, "%m")),
            as.integer(format(template, "%d")))
  })))
  dates
}

make_noleap_dates <- function(start_year, n_years) {
  # Build a template year (non-leap, e.g. 2001) and remove Feb 29 if present
  template <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  template <- template[format(template, "%m-%d") != "02-29"]  # Removes Feb 29 if present (should not be, but robust)

  dates <- as.Date(unlist(lapply(0:(n_years - 1), function(i) {
    # Use template month and day, but with incremented year
    sprintf("%04d-%02d-%02d",
            start_year + i,
            as.integer(format(template, "%m")),
            as.integer(format(template, "%d")))
  })))
  dates
}
