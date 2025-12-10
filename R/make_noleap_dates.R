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

  # Convert input
  start_date <- as.Date(start_date)

  # Build a no-leap-year template (e.g., 2001) and remove Feb 29
  template <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  template <- template[format(template, "%m-%d") != "02-29"]

  # Extract MM-DD pattern
  md <- format(template, "%m-%d")

  # Identify the month-day of the start date
  start_md <- format(start_date, "%m-%d")

  # Normalize Feb 29 ??? Feb 28
  if (start_md == "02-29") start_md <- "02-28"

  # Starting index in the no-leap template
  start_idx <- match(start_md, md)

  # Build repeating no-leap month-day sequence
  cycle_idx <- ((start_idx - 1 + 0:(n_dates - 1)) %% length(md)) + 1
  md_seq <- md[cycle_idx]

  # Compute year increments (1 year per 365 days)
  years_passed <- (0:(n_dates - 1)) %/% length(md)
  years_final <- as.integer(format(start_date, "%Y")) + years_passed

  # Construct final Date vector
  as.Date(paste0(years_final, "-", md_seq))
}
