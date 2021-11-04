

#' Calculate water years from a date series
#'
#' @param date a vector of values class of date object
#' @param water.year.first.month index of the first month in the water year,
#' e.g., if water year starts in June, it is 6.
#'
#' @return A vector of water years
#'
#'
#'
#' @export
#' @examples
#' date <- seq(as.Date("2000/1/1"), by = "month", length.out = 120)
#'
#' Water year equals to calendar year
#' getWaterYear(date, water.year.first.month = 1)
#'
#' Water year starts in October
#' getWaterYear(date, water.year.first.month = 10)
getWaterYear <- function(date = NULL, water.year.first.month = 1)  {

  w.year <- as.numeric(format(date, "%Y"))

  if(water.year.first.month != 1) {

    month.list <- c(water.year.first.month:12,1:(water.year.first.month-1))
    months_ny <- as.numeric(format(date, "%m")) > month.list[12]
    w.year[months_ny] <- w.year[months_ny] + 1
  }

  return(w.year)

}
