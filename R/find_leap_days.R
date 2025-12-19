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
