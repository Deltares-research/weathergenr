#' Potential evapotranspiration (PET) by hargreaves method
#'
#' \url{http://www.civil.uwaterloo.ca/watflood/manual/02_03_2.htm}
#' @param months To be completed...
#' @param temp To be completed...
#' @param tdiff To be completed...
#' @param lat To be completed...
#'
#' @return the output is a vector of PET values
#' @export
hargreavesPet <- function(months, temp, tdiff, lat) {

  days.j = c(15, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350)

  dr = (1 + 0.033 * cos(2 * pi/365 * days.j[months]))
  phi = pi/180 * lat
  delta = 0.409 * sin((2 * pi/365 * days.j[months]) - 1.39)
  ws = acos(-tan(phi) * tan(delta))
  Rs = ((24 * 60/pi) * 0.082 * dr * (ws * sin(phi) * sin(delta) + cos(phi) *
      cos(delta) * sin(ws))) * 0.408

  return(0.0023 * Rs * (temp + 17.8) * sqrt(tdiff))
}
