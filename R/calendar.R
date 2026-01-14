# ==============================================================================
# DATE / CALENDAR HELPERS
# ==============================================================================

#' Generate a no-leap (365-day) daily date sequence
#'
#' Constructs a continuous sequence of dates starting from a given start date,
#' using a no-leap calendar (i.e., February 29 removed). The function repeats a
#' 365-day month-day template and increments the year every 365 days.
#'
#' @param start_date Date or character coercible to Date.
#'   The first date of the sequence. If the date is February 29, it is mapped to
#'   February 28 because the no-leap template does not include Feb 29.
#' @param n_days Integer. Total number of dates to generate.
#'
#' @return A vector of class \code{Date} with length \code{n_days}, following
#'   a no-leap (365-day) calendar.
#'
#' @examples
#' generate_noleap_dates("1980-01-01", 365)
#' generate_noleap_dates(as.Date("1995-02-28"), 1000)
#'
#' @export
generate_noleap_dates <- function(start_date, n_days) {

  start_date <- as.Date(start_date)

  if (!inherits(start_date, "Date") || length(start_date) != 1L || is.na(start_date)) {
    stop("'start_date' must be a valid Date (or coercible to Date).", call. = FALSE)
  }

  if (!is.numeric(n_days) || length(n_days) != 1L || !is.finite(n_days) ||
      n_days < 1 || (n_days %% 1) != 0) {
    stop("'n_days' must be a positive integer.", call. = FALSE)
  }
  n_days <- as.integer(n_days)

  # Build no-leap month-day template (365 entries; excludes Feb 29)
  noleap_template <- seq.Date(
    as.Date("2001-01-01"),
    as.Date("2001-12-31"),
    by = "day"
  )
  noleap_template <- noleap_template[format(noleap_template, "%m-%d") != "02-29"]
  month_day <- format(noleap_template, "%m-%d")  # length 365

  # Normalize start month-day (Feb 29 -> Feb 28)
  start_month_day <- format(start_date, "%m-%d")
  if (start_month_day == "02-29") start_month_day <- "02-28"

  start_doy_noleap <- match(start_month_day, month_day)

  # Synthetic DOY sequence on no-leap calendar
  doy_seq <- ((start_doy_noleap - 1 + seq_len(n_days) - 1) %% 365) + 1

  # Detect where DOY wraps from 365 -> 1 (year increment)
  year_wrap <- c(FALSE, doy_seq[-1] < doy_seq[-length(doy_seq)])
  year_offset <- cumsum(year_wrap)

  base_year <- as.integer(format(start_date, "%Y"))
  years <- base_year + year_offset

  dates_out <- as.Date(
    paste0(years, "-", month_day[doy_seq]),
    format = "%Y-%m-%d"
  )

  dates_out
}



#' Calculate water year from a date vector
#'
#' Assigns a water year to each date based on a custom starting month.
#' A water year groups months across calendar years (e.g., Oct-Sep).
#'
#' @param date A vector of \code{Date} objects.
#' @param water_year_start_month Integer (1-12) indicating the first month of the water year.
#'   For example, \code{10} for October or \code{6} for June. Default is \code{1} (calendar year).
#'
#' @return An integer vector of the same length as \code{date}, giving the water year of each date.
#'
#' @examples
#' dates <- as.Date(c("2022-09-30", "2022-10-01", "2023-06-15"))
#' compute_water_year(dates, water_year_start_month = 10)
#'
#' @export
compute_water_year <- function(date, water_year_start_month = 1) {

  if (!inherits(date, "Date")) {
    stop("'date' must be of class Date", call. = FALSE)
  }

  if (!is.numeric(water_year_start_month) || length(water_year_start_month) != 1L ||
      !is.finite(water_year_start_month) || water_year_start_month < 1 ||
      water_year_start_month > 12 || (water_year_start_month %% 1) != 0) {
    stop("'water_year_start_month' must be an integer between 1 and 12", call. = FALSE)
  }
  water_year_start_month <- as.integer(water_year_start_month)

  cal_year  <- as.integer(format(date, "%Y"))
  cal_month <- as.integer(format(date, "%m"))

  water_year_out <- cal_year
  if (water_year_start_month > 1L) {
    idx <- cal_month >= water_year_start_month
    water_year_out[idx] <- water_year_out[idx] + 1L
  }

  as.integer(water_year_out)
}



#' Find leap-day (Feb 29) indices in a date vector
#'
#' Scans a vector of dates and returns the integer indices of all elements that
#' correspond to February 29. If the vector contains no leap days, returns \code{NULL}.
#'
#' @param date A vector coercible to \code{Date}.
#'
#' @return Integer vector of indices where dates equal February 29, or \code{NULL} if none.
#'
#' @examples
#' find_leap_day_indices(as.Date(c("1980-02-28", "1980-02-29", "1981-01-01")))
#' find_leap_day_indices(as.Date(c("2001-01-01", "2001-12-31")))
#'
#' @export
find_leap_day_indices <- function(date) {

  date <- as.Date(date)
  if (anyNA(date)) {
    stop("'date' must be coercible to Date with no NA values.", call. = FALSE)
  }

  leap_idx <- which(format(date, "%m-%d") == "02-29")

  if (length(leap_idx) == 0L) NULL else leap_idx
}
