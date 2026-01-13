#' @title Estimate Monthly Potential Evapotranspiration (Hargreaves Method)
#'
#' @description
#' Computes monthly mean potential evapotranspiration (PET) using the
#' Hargreaves temperature-based method. PET is estimated from mean air
#' temperature, diurnal temperature range, and extraterrestrial radiation
#' derived from latitude and time of year.
#'
#' The implementation follows the FAO-56 formulation and returns PET
#' in units of millimeters per day (mm/day).
#'
#' @param months Integer vector (1-12). Calendar months corresponding
#'   to each PET estimate.
#' @param temp Numeric vector. Mean monthly air temperature (degrees C).
#' @param tdiff Numeric vector. Monthly diurnal temperature range
#'   (Tmax ??? Tmin, degrees C). Must be non-negative.
#' @param lat Numeric scalar. Latitude in decimal degrees
#'   (positive north, negative south).
#'
#' @return
#' Numeric vector of monthly potential evapotranspiration values
#' in mm/day, with length equal to `length(months)`.
#'
#' @details
#' Extraterrestrial radiation is computed using a fixed representative
#' day of year for the middle of each month. Radiation is converted from
#' MJ m^2/day to equivalent evaporation depth using the FAO conversion
#' factor of 0.408.
#'
#' The Hargreaves equation is:
#' \deqn{
#' PET = 0.0023 * Ra * (T + 17.8) * sqrt(Tmax - Tmin)
#' }
#' where Ra is extraterrestrial radiation.
#'
#' @references
#' Hargreaves, G.H. and Samani, Z.A. (1985).
#' Reference crop evapotranspiration from temperature.
#' Applied Engineering in Agriculture, 1(2), 96-99.
#'
#' FAO (1998). Crop Evapotranspiration (FAO Irrigation and Drainage Paper 56).
#'
#' @export
pet_hargreaves <- function(months, temp, tdiff, lat) {

  ## ----------------------------
  ## Input validation
  ## ----------------------------
  if (!is.numeric(months)) {
    stop("'months' must be numeric")
  }

  months <- as.integer(months)
  if (anyNA(months) || any(months < 1 | months > 12)) {
    stop("'months' must contain integers between 1 and 12")
  }

  n <- length(months)
  if (length(temp) != n || length(tdiff) != n) {
    stop("'months', 'temp', and 'tdiff' must have the same length")
  }

  if (!is.numeric(temp) || any(!is.finite(temp))) {
    stop("'temp' must contain finite numeric values")
  }

  if (!is.numeric(tdiff) || any(tdiff < 0, na.rm = TRUE)) {
    stop("'tdiff' must be non-negative")
  }

  if (!is.numeric(lat) || length(lat) != 1L || !is.finite(lat) ||
      lat < -90 || lat > 90) {
    stop("'lat' must be a finite numeric value between -90 and 90")
  }

  ## ----------------------------
  ## Day-of-year approximation
  ## ----------------------------
  doy_mid_month <- c(15, 46, 75, 106, 136, 167,
                     197, 228, 259, 289, 320, 350)
  doy <- doy_mid_month[months]

  ## ----------------------------
  ## Solar geometry
  ## ----------------------------
  lat_rad <- lat * pi / 180

  inv_rel_dist <- 1 + 0.033 * cos(2 * pi * doy / 365)
  solar_decl   <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  sunset_angle <- acos(-tan(lat_rad) * tan(solar_decl))

  ## Extraterrestrial radiation (MJ m-2 day-1)
  ra <- (24 * 60 / pi) * 0.082 * inv_rel_dist *
    (sunset_angle * sin(lat_rad) * sin(solar_decl) +
       cos(lat_rad) * cos(solar_decl) * sin(sunset_angle))

  ## Convert radiation to equivalent evaporation (mm/day)
  ra_mm <- 0.408 * ra

  ## ----------------------------
  ## Hargreaves PET
  ## ----------------------------
  0.0023 * ra_mm * (temp + 17.8) * sqrt(tdiff)

}

