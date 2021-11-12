

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
getWaterYear <- function(date = NULL, water.year.first.month = 1)  {

  w.year <- as.numeric(format(date, "%Y"))

  if(water.year.first.month != 1) {

    month.list <- c(water.year.first.month:12,1:(water.year.first.month-1))
    months_ny <- as.numeric(format(date, "%m")) > month.list[12]
    w.year[months_ny] <- w.year[months_ny] + 1
  }

  return(w.year)

}
