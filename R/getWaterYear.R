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
#' getWaterYear(dates, water.year.first.month = 10)
#' # [1] 2022 2023 2023
#'
#' # Water year starting in June (e.g., for monsoon-driven systems)
#' dates <- as.Date(c("2023-05-31", "2023-06-01", "2023-12-01"))
#' getWaterYear(dates, water.year.first.month = 6)
#' # [1] 2023 2024 2024
#'
#' # Default: calendar year (Jan-Dec)
#' dates <- as.Date(c("2024-01-01", "2024-12-31"))
#' getWaterYear(dates)
#' # [1] 2024 2024
#'
#' @export
getWaterYear <- function(date, water.year.first.month = 1) {
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
