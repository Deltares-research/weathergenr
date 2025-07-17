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

  # Make a sequence of years
  years <- seq(start_year, by = 1, length.out = n_years)

  # Extract month and day for template year
  month_day <- format(template, "%m-%d")

  # Efficiently build all dates for each year
  # Outer product: years (rows) x days (cols)
  all_dates <- as.Date(paste0(
    rep(years, each = length(month_day)), "-",
    rep(month_day, times = n_years)
  ))

  return(all_dates)
}
