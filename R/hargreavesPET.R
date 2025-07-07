#' @title Estimate Potential Evapotranspiration Using the Hargreaves Method
#'
#' @description
#' Computes monthly potential evapotranspiration (PET) using the Hargreaves method,
#' which estimates PET from temperature and extraterrestrial radiation.
#'
#' Reference: \url{http://www.civil.uwaterloo.ca/watflood/manual/02_03_2.htm}
#'
#' @param months An integer vector of month indices (1-12).
#' @param temp A numeric vector of mean monthly temperatures in degrees Celsius.
#' @param tdiff A numeric vector of monthly temperature differences (Tmax - Tmin) in degrees Celsius.
#' @param lat A single numeric value representing latitude in decimal degrees.
#'
#' @return A numeric vector of PET values (in mm/day), one for each month.
#' @examples
#' # Example: Estimate PET for three months in Amsterdam (lat 52.37)
#' months <- c(1, 6, 12)                         # January, June, December
#' mean_temp <- c(3.2, 16.1, 4.0)                # Example monthly mean temperatures (DegC)
#' temp_range <- c(6.0, 8.5, 5.2)                # Example monthly Tmax-Tmin (DegC)
#' latitude <- 52.37                             # Amsterdam
#' pet_values <- hargreavesPET(months, mean_temp, temp_range, latitude)
#' print(pet_values)
#'
#' # Estimate PET for a full year at 40°N
#' months <- 1:12
#' mean_temp <- c(-1.1, 0.5, 5.6, 11.2, 16.0, 20.2, 22.7, 22.3, 18.2, 12.0, 5.4, 0.3)
#' temp_range <- c(5, 6, 7, 8, 9, 10, 10, 9, 8, 7, 6, 5)
#' latitude <- 40
#' pet_year <- hargreavesPET(months, mean_temp, temp_range, latitude)
#' plot(months, pet_year, type = "o", xlab = "Month", ylab = "PET (mm/day)")
#'
#' @export
hargreavesPET <- function(months, temp, tdiff, lat) {

  # Validate inputs
  if (any(months < 1 | months > 12)) {
    stop("All 'months' values must be between 1 and 12.")
  }
  if (length(lat) != 1 || !is.numeric(lat) || lat < -90 || lat > 90) {
    stop("'lat' must be a numeric value between -90 and 90.")
  }

  # Day of year at the middle of each month (approximation)
  days.j <- c(15, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350)
  j <- days.j[months]

  # Convert latitude to radians
  phi <- pi / 180 * lat

  # Calculate solar-related components
  dr <- 1 + 0.033 * cos(2 * pi / 365 * j)         # Inverse relative distance Earth-Sun
  delta <- 0.409 * sin(2 * pi / 365 * j - 1.39)   # Solar declination
  ws <- acos(-tan(phi) * tan(delta))              # Sunset hour angle

  # Extraterrestrial radiation (Ra) in MJ/m/day, converted to mm/day with 0.408
  Ra <- (24 * 60 / pi) * 0.082 * dr * (ws * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(ws))
  Rs <- 0.408 * Ra

  # Hargreaves PET formula
  PET <- 0.0023 * Rs * (temp + 17.8) * sqrt(tdiff)

  return(PET)
}
