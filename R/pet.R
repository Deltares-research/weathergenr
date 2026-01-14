#' @title Calculate Monthly Potential Evapotranspiration (PET)
#'
#' @description
#' Calculates monthly potential evapotranspiration (PET) using a selected method.
#' Currently supports the Hargreaves temperature-based method.
#'
#' @param month Integer vector (1-12). Calendar month for each PET estimate.
#' @param temp Numeric vector. Mean monthly air temperature (degC).
#' @param temp_range Numeric vector. Monthly diurnal temperature range (Tmax - Tmin, degC).
#'   Must be non-negative.
#' @param lat_deg Numeric scalar. Latitude in decimal degrees (positive north).
#' @param method Character scalar. PET method. Currently \code{"hargreaves"}.
#'
#' @return Numeric vector of monthly PET values in mm/day (length \code{length(month)}).
#'
#' @export
calculate_monthly_pet <- function(
    month,
    temp,
    temp_range,
    lat_deg,
    method = "hargreaves"
) {

  # ---------------------------------------------------------------------------
  # Method dispatch
  # ---------------------------------------------------------------------------
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("'method' must be a single non-missing character value", call. = FALSE)
  }
  method <- tolower(method)

  if (!method %in% c("hargreaves")) {
    stop(
      "'method' must be one of: ",
      paste(c("hargreaves"), collapse = ", "),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Shared input validation (applies to all methods)
  # ---------------------------------------------------------------------------
  if (!is.numeric(month)) stop("'month' must be numeric", call. = FALSE)
  month <- as.integer(month)
  if (anyNA(month) || any(month < 1L | month > 12L)) {
    stop("'month' must contain integers between 1 and 12", call. = FALSE)
  }

  n <- length(month)
  if (length(temp) != n || length(temp_range) != n) {
    stop("'month', 'temp', and 'temp_range' must have the same length", call. = FALSE)
  }

  if (!is.numeric(temp) || any(!is.finite(temp))) {
    stop("'temp' must contain finite numeric values", call. = FALSE)
  }

  if (!is.numeric(temp_range) || any(!is.finite(temp_range)) || any(temp_range < 0)) {
    stop("'temp_range' must contain finite non-negative numeric values", call. = FALSE)
  }

  if (!is.numeric(lat_deg) || length(lat_deg) != 1L || !is.finite(lat_deg) ||
      lat_deg < -90 || lat_deg > 90) {
    stop("'lat_deg' must be a finite numeric value between -90 and 90", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Dispatch to implementation
  # ---------------------------------------------------------------------------
  if (method == "hargreaves") {
    return(.pet_monthly_hargreaves(
      month = month,
      temp = temp,
      temp_range = temp_range,
      lat_deg = lat_deg
    ))
  }

  stop("Internal error: unsupported 'method' dispatch", call. = FALSE)
}


#' Monthly PET via Hargreaves (FAO-56)
#' @keywords internal
.pet_monthly_hargreaves <- function(month, temp, temp_range, lat_deg) {

  doy_mid <- c(15, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350)
  doy <- doy_mid[month]

  lat_rad <- lat_deg * pi / 180

  dr <- 1 + 0.033 * cos(2 * pi * doy / 365)
  delta <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  omega_s <- acos(-tan(lat_rad) * tan(delta))

  ra_mj_m2_day <- (24 * 60 / pi) * 0.082 * dr *
    (omega_s * sin(lat_rad) * sin(delta) +
       cos(lat_rad) * cos(delta) * sin(omega_s))

  ra_mm_day <- 0.408 * ra_mj_m2_day

  0.0023 * ra_mm_day * (temp + 17.8) * sqrt(temp_range)
}
