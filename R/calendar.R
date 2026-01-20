# ==============================================================================
# DATE / CALENDAR HELPERS
# ==============================================================================

#' Generate a no-leap (365-day) daily date sequence
#'
#' Constructs a continuous sequence of dates starting from a given start date,
#' using a no-leap calendar (i.e., February 29 removed). The function repeats a
#' 365-day month-day template and increments the year every 365 days.
#'
#' @param start_date Date or character coercible to Date.
#'   The first date of the sequence. If the date is February 29, it is mapped to
#'   February 28 because the no-leap template does not include Feb 29.
#' @param n_days Integer. Total number of dates to generate.
#'
#' @return A vector of class \code{Date} with length \code{n_days}, following
#'   a no-leap (365-day) calendar.
#'
#' @examples
#' generate_noleap_dates("1980-01-01", 365)
#' generate_noleap_dates(as.Date("1995-02-28"), 1000)
#'
#' @export
generate_noleap_dates <- function(start_date, n_days) {

  start_date <- .safe_as_date(start_date)

  if (!inherits(start_date, "Date") || length(start_date) != 1L || is.na(start_date)) {
    stop("'start_date' must be a valid Date (or coercible to Date).", call. = FALSE)
  }

  if (!is.numeric(n_days) || length(n_days) != 1L || !is.finite(n_days) ||
      n_days < 1 || (n_days %% 1) != 0) {
    stop("'n_days' must be a positive integer.", call. = FALSE)
  }
  n_days <- as.integer(n_days)

  # Build no-leap month-day template (365 entries; excludes Feb 29)
  noleap_template <- seq.Date(
    as.Date("2001-01-01"),
    as.Date("2001-12-31"),
    by = "day"
  )
  noleap_template <- noleap_template[format(noleap_template, "%m-%d") != "02-29"]
  month_day <- format(noleap_template, "%m-%d")  # length 365

  # Normalize start month-day (Feb 29 -> Feb 28)
  start_month_day <- format(start_date, "%m-%d")
  if (start_month_day == "02-29") start_month_day <- "02-28"

  start_doy_noleap <- match(start_month_day, month_day)

  # Synthetic DOY sequence on no-leap calendar
  doy_seq <- ((start_doy_noleap - 1 + seq_len(n_days) - 1) %% 365) + 1

  # Detect where DOY wraps from 365 -> 1 (year increment)
  year_wrap <- c(FALSE, doy_seq[-1] < doy_seq[-length(doy_seq)])
  year_offset <- cumsum(year_wrap)

  base_year <- as.integer(format(start_date, "%Y"))
  years <- base_year + year_offset

  dates_out <- as.Date(
    paste0(years, "-", month_day[doy_seq]),
    format = "%Y-%m-%d"
  )

  dates_out
}



#' Calculate water year from a date vector
#'
#' Assigns a water year to each date based on a custom starting month.
#' A water year groups months across calendar years (e.g., Oct-Sep).
#'
#' @param date A vector of \code{Date} objects.
#' @param water_year_start_month Integer (1-12) indicating the first month of the water year.
#'   For example, \code{10} for October or \code{6} for June. Default is \code{1} (calendar year).
#'
#' @return An integer vector of the same length as \code{date}, giving the water year of each date.
#'
#' @examples
#' dates <- as.Date(c("2022-09-30", "2022-10-01", "2023-06-15"))
#' compute_water_year(dates, water_year_start_month = 10)
#'
#' @export
compute_water_year <- function(date, water_year_start_month = 1) {

  if (!inherits(date, "Date")) {
    stop("'date' must be of class Date", call. = FALSE)
  }

  if (!is.numeric(water_year_start_month) || length(water_year_start_month) != 1L ||
      !is.finite(water_year_start_month) || water_year_start_month < 1 ||
      water_year_start_month > 12 || (water_year_start_month %% 1) != 0) {
    stop("'water_year_start_month' must be an integer between 1 and 12", call. = FALSE)
  }
  water_year_start_month <- as.integer(water_year_start_month)

  cal_year  <- as.integer(format(date, "%Y"))
  cal_month <- as.integer(format(date, "%m"))

  water_year_out <- cal_year
  if (water_year_start_month > 1L) {
    idx <- cal_month >= water_year_start_month
    water_year_out[idx] <- water_year_out[idx] + 1L
  }

  as.integer(water_year_out)
}



#' Find leap-day (Feb 29) indices in a date vector
#'
#' Scans a vector of dates and returns the integer indices of all elements that
#' correspond to February 29. If the vector contains no leap days, returns \code{NULL}.
#'
#' @param date A vector coercible to \code{Date}.
#'
#' @return Integer vector of indices where dates equal February 29, or \code{NULL} if none.
#'
#' @examples
#' find_leap_day_indices(as.Date(c("1980-02-28", "1980-02-29", "1981-01-01")))
#' find_leap_day_indices(as.Date(c("2001-01-01", "2001-12-31")))
#'
#' @export
find_leap_day_indices <- function(date) {

  date <- .safe_as_date(date)

  if (anyNA(date)) {
    stop("'date' must be coercible to Date with no NA values.", call. = FALSE)
  }

  leap_idx <- which(format(date, "%m-%d") == "02-29")

  if (length(leap_idx) == 0L) NULL else leap_idx
}


# ==============================================================================
# INTERNAL HELPER FUNCTIONS FOR generate_weather()
# ==============================================================================

#' Normalize Calendar to 365 Days
#'
#' @description
#' Removes February 29 (leap days) from observation dates and corresponding
#' rows from all grid cell data frames to enforce a 365-day calendar.
#'
#' @param obs_dates Date vector of observation dates.
#' @param obs_data Named list of data frames (one per grid cell).
#' @param verbose Logical. If TRUE, logs the number of dropped days.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{dates}{Date vector with leap days removed.}
#'   \item{data}{List of data frames with corresponding rows removed.}
#' }
#'
#' @export
normalize_calendar <- function(obs_dates, obs_data, verbose = FALSE) {

  leap_idx <- find_leap_day_indices(obs_dates)

  if (is.null(leap_idx) || length(leap_idx) == 0L) {
    return(list(dates = obs_dates, data = obs_data))
  }

  n_orig <- length(obs_dates)
  dates_out <- obs_dates[-leap_idx]

  data_out <- lapply(obs_data, function(df) {
    if (!is.data.frame(df)) {
      stop("Each element of obs_data must be a data.frame.", call. = FALSE)
    }
    if (nrow(df) != n_orig) {
      stop("Each obs_data element must have nrow equal to length(obs_dates).",
           call. = FALSE)
    }
    df[-leap_idx, , drop = FALSE]
  })

  .log(
    "Dropped {length(leap_idx)} leap days (365-day calendar enforced)",
    tag = "INIT",
    verbose = verbose
  )

  list(dates = dates_out, data = data_out)
}


#' Build Historical Date Table
#'
#' @description
#' Constructs a tibble of historical dates with calendar and water-year indices,
#' then restricts to complete years (365 days each) and returns the longest
#' contiguous block of complete years.
#'
#' @param obs_dates Date vector of observation dates (leap days already removed).
#' @param year_start_month Integer 1-12. First month of the simulation year.
#' @param verbose Logical. If TRUE, logs summary information.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{dates_df}{Tibble with columns: dateo, year, wyear, month, day, date.}
#'   \item{wyear_idx}{Integer vector of row indices in obs_dates corresponding
#'     to complete years.}
#'   \item{complete_wyears}{Integer vector of complete water years.}
#' }
#'
#' @export
build_historical_dates <- function(obs_dates, year_start_month, verbose = FALSE) {

  his_wyear <- compute_water_year(obs_dates, year_start_month)

  dates_df <- tibble::tibble(
    dateo = obs_dates,
    year  = as.integer(format(obs_dates, "%Y")),
    wyear = his_wyear,
    month = as.integer(format(obs_dates, "%m")),
    day   = as.integer(format(obs_dates, "%d"))
  )

  # Compute internal date representation
  if (year_start_month == 1L) {
    dates_df$date <- dates_df$dateo
  } else {
    dates_df$date <- as.Date(
      sprintf("%04d-%02d-%02d", dates_df$wyear, dates_df$month, dates_df$day)
    )
  }

  # Validation
  if (anyNA(dates_df$date)) {
    stop("Internal error: dates_df$date contains NA.", call. = FALSE)
  }
  if (anyDuplicated(dates_df$date)) {
    stop("Internal error: dates_df$date contains duplicates.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Identify complete years (365 days) and keep longest contiguous block
  # ---------------------------------------------------------------------------
  wyear_counts <- as.data.frame(table(dates_df$wyear), stringsAsFactors = FALSE)
  names(wyear_counts) <- c("wyear", "n_days")
  wyear_counts$wyear <- as.integer(wyear_counts$wyear)

  full_wyears <- wyear_counts$wyear[wyear_counts$n_days == 365L]

  if (length(full_wyears) == 0L) {
    stop(
      "No complete years found (need 365 days per year after leap day removal).",
      call. = FALSE
    )
  }

  # Find longest contiguous block
  full_wyears <- sort(full_wyears)
  runs <- split(full_wyears, cumsum(c(1L, diff(full_wyears) != 1L)))
  full_wyears <- runs[[which.max(lengths(runs))]]

  # Filter to complete years
  dates_df <- dplyr::filter(dates_df, wyear %in% full_wyears)

  # Build index mapping back to original obs_dates
  wyear_idx <- match(dates_df$dateo, obs_dates)

  if (anyNA(wyear_idx)) {
    stop("Internal error: dates_df$dateo did not match obs_dates.", call. = FALSE)
  }

  .log(
    "Complete years: {length(full_wyears)} ({min(full_wyears)}-{max(full_wyears)})",
    tag = "INIT",
    verbose = verbose
  )

  list(
    dates_df       = dates_df,
    wyear_idx      = wyear_idx,
    complete_wyears = full_wyears
  )
}



# ==============================================================================
# INTERNAL: Safe Date coercion (turns as.Date() errors into NA)
# ==============================================================================

.safe_as_date <- function(x) {
  if (inherits(x, "Date")) return(x)

  tryCatch(
    as.Date(x),
    error = function(e) {
      # Preserve length for vector inputs
      rep(as.Date(NA), length(x))
    }
  )
}

