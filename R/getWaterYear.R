

#' Water year function
#'
#' @param month.list placeholder
#' @param date placeholder
#'
#' @return
#' @export
getWaterYear <- function(date = NULL, month.list = 1:12)  {

  w.year <- as.numeric(format(date, "%Y"))

  if (all(month.list != 1:12)) {
    months_ny <- as.numeric(format(date, "%m")) > month.list[12]
    w.year[months_ny] <- w.year[months_ny] + 1
  }

  return(w.year)
}
