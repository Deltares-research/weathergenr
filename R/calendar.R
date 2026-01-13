#' Generate a no-leap (365-day) daily date sequence
#'
#' Constructs a continuous sequence of dates starting from a given start date,
#' using a no-leap calendar (i.e., February 29 removed). The function repeats a
#' 365-day month-day template and increments the year every 365 days.
#'
#' @param start_date Date or character coercible to Date.
#'   The first date of the sequence. If the date is February 29, it is mapped to
#'   February 28 because the no-leap template does not include Feb 29.
#' @param n_dates Integer. Total number of dates to generate.
#'
#' @return A vector of class \code{Date} with length \code{n_dates}, following
#'   a no-leap (365-day) calendar.
#'
#' @examples
#' make_noleap_dates("1980-01-01", 365)
#' make_noleap_dates(as.Date("1995-02-28"), 1000)
#'
#' @export
make_noleap_dates <- function(start_date, n_dates) {

  start_date <- as.Date(start_date)

  # Build no-leap template of 365 day-of-year ??? MM-DD
  tpl <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  tpl <- tpl[format(tpl, "%m-%d") != "02-29"]
  md <- format(tpl, "%m-%d")               # vector of length 365

  # Identify start MD (normalize Feb 29 ??? Feb 28)
  start_md <- format(start_date, "%m-%d")
  if (start_md == "02-29") start_md <- "02-28"

  start_doy <- match(start_md, md)         # synthetic DOY

  # Synthetic DOY sequence
  seq_doy <- ((start_doy - 1 + seq_len(n_dates) - 1) %% 365) + 1

  # Detect where DOY wraps from 365 ??? 1
  wraps <- c(FALSE, seq_doy[-1] < seq_doy[-length(seq_doy)])

  # Cumulative year increments
  year_inc <- cumsum(wraps)

  # Base year
  years <- as.integer(format(start_date, "%Y")) + year_inc

  # Construct dates
  out <- as.Date(paste0(years, "-", md[seq_doy]), format = "%Y-%m-%d")

  out
}



#' @title Calculate Water Year from Date Vector
#'
#' @description
#' Assigns a water year to each date in a vector based on a custom starting month.
#' Common in hydrology, a water year groups months across calendar years (e.g., Oct-Sep).
#'
#' @param date A vector of `Date` objects.
#' @param water.year.first.month Integer (1-12) indicating the first month of the water year.
#' For example, `10` for October or `6` for June.
#'
#' @return An integer vector of the same length as `date`, giving the water year of each date.
#'
#' #' @examples
#' # Water year starting in October (common in US hydrology)
#' dates <- as.Date(c("2022-09-30", "2022-10-01", "2023-06-15"))
#' get_water_year(dates, water.year.first.month = 10)
#' # [1] 2022 2023 2023
#'
#' # Water year starting in June (e.g., for monsoon-driven systems)
#' dates <- as.Date(c("2023-05-31", "2023-06-01", "2023-12-01"))
#' get_water_year(dates, water.year.first.month = 6)
#' # [1] 2023 2024 2024
#'
#' # Default: calendar year (Jan-Dec)
#' dates <- as.Date(c("2024-01-01", "2024-12-31"))
#' get_water_year(dates)
#' # [1] 2024 2024#'
#' @export
get_water_year <- function(date, water.year.first.month = 1) {
  # Input validation
  if (!inherits(date, "Date")) {
    stop("'date' must be of class Date")
  }

  if (!is.numeric(water.year.first.month) || water.year.first.month < 1 || water.year.first.month > 12) {
    stop("'water.year.first.month' must be an integer between 1 and 12")
  }

  # Get calendar year and month
  year <- as.integer(format(date, "%Y"))
  month <- as.integer(format(date, "%m"))

  # Adjust year for months >= start of water year
  water_year <- year
  if (water.year.first.month > 1) {
    water_year[month >= water.year.first.month] <- water_year[month >= water.year.first.month] + 1
  }
  return(as.integer(water_year))
}

#' Find leap-day positions (Feb 29) in a date vector
#'
#' Scans a vector of dates and returns the integer indices of all elements that
#' correspond to February 29. If the vector contains no leap days, the function
#' returns \code{NULL}.
#'
#' @param dates A vector coercible to \code{Date}.
#'
#' @return Integer vector of indices where dates equal February 29, or
#'   \code{NULL} if no leap days are present.
#'
#' @examples
#' find_leap_days(as.Date(c("1980-02-28", "1980-02-29", "1981-01-01")))
#' # Returns: 2
#'
#' find_leap_days(as.Date(c("2001-01-01", "2001-12-31")))
#' # Returns: NULL
#'
#' @export
find_leap_days <- function(dates) {
  dates <- as.Date(dates)
  idx <- which(format(dates, "%m-%d") == "02-29")

  if (length(idx) == 0) {
    return(NULL)
  } else {
    return(idx)
  }
}
