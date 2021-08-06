#' Potential evapotranspiration (PET) by hargreaves method
#'
#' \url{http://www.civil.uwaterloo.ca/watflood/manual/02_03_2.htm}
#'
#' @param months a time-series date object
#' @param temp a vector of average temperature values (°C)
#' @param tdif a vector of differences computed from maximum and minimum
#' temperatures
#' @param lat latitude information (negative values for southern hemisphere)

#' @return the output is a vector of PET values
#' @export
hargreavesPET <- function(months, temp, tdiff, lat) {

  days.m = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  days.j = c(15, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350)

  mday <- days.m[months]

  dr = (1 + 0.033 * cos(2 * pi/365 * days.j[months]))
  phi = pi/180 * lat
  delta = 0.409 * sin((2 * pi/365 * days.j[months]) - 1.39)
  ws = acos(-tan(phi) * tan(delta))
  Rs = ((24 * 60/pi) * 0.082 * dr * (ws * sin(phi) * sin(delta) + cos(phi) *
      cos(delta) * sin(ws))) * 0.408 * mday

  return(0.0023 * Rs * (temp + 17.8) * sqrt(tdiff))
}
