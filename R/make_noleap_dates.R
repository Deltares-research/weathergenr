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

