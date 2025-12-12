#' Calculate Monthly Markov Transition Probabilities for a 3-State Chain
#'
#' @description
#' Computes month-specific transition probabilities for a three-state
#' precipitation occurrence Markov chain (dry, wet, very wet) using
#' lagged daily precipitation and month-dependent thresholds.
#'
#' The function supports both calendar-year and water-year simulations.
#' The regime is inferred automatically from `month_list`:
#' if `month_list[1] == 1`, calendar-year logic is used; otherwise,
#' water-year logic is assumed.
#'
#' Transition probabilities are estimated from observed lag-1 to lag-0
#' precipitation pairs, optionally filtered to enforce same-year lags
#' (recommended for calendar-year mode). Probabilities are then adjusted
#' for dry- and wet-spell persistence and assigned to the appropriate
#' simulated days of the target simulation year.
#'
#' @param PRCP_LAG0 Numeric vector. Precipitation at lag 0 (current day)
#'   for observed transitions.
#' @param PRCP_LAG1 Numeric vector. Precipitation at lag 1 (previous day)
#'   for observed transitions.
#' @param MONTH_LAG0 Integer vector (1-12). Calendar month corresponding
#'   to `PRCP_LAG0`.
#' @param MONTH_LAG1 Integer vector (1-12). Calendar month corresponding
#'   to `PRCP_LAG1`.
#' @param YEAR_LAG0 Optional integer vector. Year or water-year
#'   corresponding to `PRCP_LAG0`. Used to enforce same-year lag filtering.
#' @param YEAR_LAG1 Optional integer vector. Year or water-year
#'   corresponding to `PRCP_LAG1`. Used to enforce same-year lag filtering.
#' @param wet_threshold Numeric vector of length 12. Monthly precipitation
#'   thresholds separating dry and wet states, aligned to `month_list`.
#' @param extreme_threshold Numeric vector of length 12. Monthly
#'   precipitation thresholds separating wet and very wet states, aligned
#'   to `month_list`.
#' @param month_list Integer vector of length 12. Defines the simulated
#'   year ordering of months (for example 1:12 for calendar year, or
#'   c(10:12, 1:9) for an October-start water year).
#' @param MONTH_SIM Integer vector. Simulation calendar months for each
#'   simulated day.
#' @param WATER_YEAR_SIM Integer vector. Simulation year (calendar or
#'   water year, depending on regime) for each simulated day.
#' @param y Integer. Index of the current simulated year (1-based).
#' @param START_YEAR_SIM Integer. Starting year of the simulation.
#'   Interpreted as a calendar year or water year depending on regime.
#' @param dry.spell.change Numeric vector of length 12. Monthly adjustment
#'   factors controlling dry-spell persistence.
#' @param wet.spell.change Numeric vector of length 12. Monthly adjustment
#'   factors controlling wet-spell persistence.
#' @param SIM_LENGTH Integer. Total number of simulated days. Determines
#'   the length of the returned probability vectors.
#'
#' @return
#' A named list of nine numeric vectors, each of length `SIM_LENGTH`:
#' \itemize{
#'   \item p00_final: P(dry to dry)
#'   \item p01_final: P(dry to wet)
#'   \item p02_final: P(dry to very wet)
#'   \item p10_final: P(wet to dry)
#'   \item p11_final: P(wet to wet)
#'   \item p12_final: P(wet to very wet)
#'   \item p20_final: P(very wet to dry)
#'   \item p21_final: P(very wet to wet)
#'   \item p22_final: P(very wet to very wet)
#' }
#'
#' @details
#' States are defined as:
#' \describe{
#'   \item{Dry}{Precipitation less than or equal to the monthly wet threshold}
#'   \item{Wet}{Precipitation greater than the wet threshold and less than
#'              or equal to the extreme threshold}
#'   \item{Very wet}{Precipitation greater than the extreme threshold}
#' }
#'
#' If `YEAR_LAG0` and `YEAR_LAG1` are provided, transitions are filtered to
#' same-year pairs. If they are not provided and calendar-year mode is
#' active, a weaker month-based guard is applied to avoid December-January
#' cross-year contamination.
#'
#' The estimated probabilities are assigned to all simulated days that fall
#' within the target simulation year and corresponding month.
#'
#' @export
calculateMarkovProbs <- function(
    PRCP_LAG0,
    PRCP_LAG1,
    MONTH_LAG0,
    MONTH_LAG1,
    YEAR_LAG0 = NULL,
    YEAR_LAG1 = NULL,
    wet_threshold,
    extreme_threshold,
    month_list,
    MONTH_SIM,
    WATER_YEAR_SIM,
    y,
    START_YEAR_SIM,
    dry.spell.change,
    wet.spell.change,
    SIM_LENGTH
) {

  # -------------------------------------------------
  # Calendar regime derived from month_list
  # month_list[1] == 1  -> calendar-year mode
  # month_list[1] != 1  -> water-year mode
  # -------------------------------------------------
  water.year <- (month_list[1] != 1)

  # -------------------------------------------------
  # Target modeled year in the simulation time axis
  # water-year mode: first WY is START_YEAR_SIM + 1
  # calendar-year mode: first year is START_YEAR_SIM
  # -------------------------------------------------
  target_year <- if (water.year) (START_YEAR_SIM + y) else (START_YEAR_SIM + y - 1)

  # -------------------------------------------------
  # Optional: strong lag guard using modeled year labels
  # -------------------------------------------------
  if (!is.null(YEAR_LAG0) && !is.null(YEAR_LAG1)) {
    same_year <- (YEAR_LAG0 == YEAR_LAG1)
    PRCP_LAG0  <- PRCP_LAG0[same_year]
    PRCP_LAG1  <- PRCP_LAG1[same_year]
    MONTH_LAG0 <- MONTH_LAG0[same_year]
    MONTH_LAG1 <- MONTH_LAG1[same_year]
  } else {
    # Fallback guard if YEAR lags are not supplied (weaker)
    if (!water.year) {
      keep <- (MONTH_LAG0 >= MONTH_LAG1)
      PRCP_LAG0  <- PRCP_LAG0[keep]
      PRCP_LAG1  <- PRCP_LAG1[keep]
      MONTH_LAG0 <- MONTH_LAG0[keep]
      MONTH_LAG1 <- MONTH_LAG1[keep]
    }
  }

  p00_final <- rep(NA_real_, SIM_LENGTH)
  p01_final <- rep(NA_real_, SIM_LENGTH)
  p02_final <- rep(NA_real_, SIM_LENGTH)

  p10_final <- rep(NA_real_, SIM_LENGTH)
  p11_final <- rep(NA_real_, SIM_LENGTH)
  p12_final <- rep(NA_real_, SIM_LENGTH)

  p20_final <- rep(NA_real_, SIM_LENGTH)
  p21_final <- rep(NA_real_, SIM_LENGTH)
  p22_final <- rep(NA_real_, SIM_LENGTH)

  # -------------------------------------------------
  # Estimate month-specific transition probabilities
  # and assign them to the simulation year "target_year"
  # -------------------------------------------------
  for (m in seq_along(month_list)) {

    mm <- month_list[m]

    r <- which(MONTH_SIM == mm & WATER_YEAR_SIM == target_year)
    if (!length(r)) next

    x <- which(MONTH_LAG1 == mm)
    if (!length(x)) {
      # No information for this month, set to 0-prob transitions
      p00_final[r] <- 1; p01_final[r] <- 0; p02_final[r] <- 0
      p10_final[r] <- 1; p11_final[r] <- 0; p12_final[r] <- 0
      p20_final[r] <- 1; p21_final[r] <- 0; p22_final[r] <- 0
      next
    }

    wet_thr <- wet_threshold[m]
    ext_thr <- extreme_threshold[m]

    # Protect against non-finite thresholds
    if (!is.finite(wet_thr)) wet_thr <- stats::quantile(PRCP_LAG1, 0.2, na.rm = TRUE)
    if (!is.finite(ext_thr)) {
      pos <- PRCP_LAG1[PRCP_LAG1 > 0]
      if (!length(pos)) pos <- PRCP_LAG1
      ext_thr <- stats::quantile(pos, 0.8, na.rm = TRUE)
    }

    # Define states for lag1 and lag0
    lag1_dry  <- (PRCP_LAG1[x] <= wet_thr)
    lag1_wet  <- (PRCP_LAG1[x] >  wet_thr & PRCP_LAG1[x] <= ext_thr)
    lag1_vwet <- (PRCP_LAG1[x] >  ext_thr)

    lag0_dry  <- (PRCP_LAG0[x] <= wet_thr)
    lag0_wet  <- (PRCP_LAG0[x] >  wet_thr & PRCP_LAG0[x] <= ext_thr)
    lag0_vwet <- (PRCP_LAG0[x] >  ext_thr)

    denom_dry  <- sum(lag1_dry)
    denom_wet  <- sum(lag1_wet)
    denom_vwet <- sum(lag1_vwet)

    p00 <- if (denom_dry  > 0) sum(lag1_dry  & lag0_dry)  / denom_dry  else 0
    p01 <- if (denom_dry  > 0) sum(lag1_dry  & lag0_wet)  / denom_dry  else 0
    p02 <- if (denom_dry  > 0) sum(lag1_dry  & lag0_vwet) / denom_dry  else 0

    p10 <- if (denom_wet  > 0) sum(lag1_wet  & lag0_dry)  / denom_wet  else 0
    p11 <- if (denom_wet  > 0) sum(lag1_wet  & lag0_wet)  / denom_wet  else 0
    p12 <- if (denom_wet  > 0) sum(lag1_wet  & lag0_vwet) / denom_wet  else 0

    p20 <- if (denom_vwet > 0) sum(lag1_vwet & lag0_dry)  / denom_vwet else 0
    p21 <- if (denom_vwet > 0) sum(lag1_vwet & lag0_wet)  / denom_vwet else 0
    p22 <- if (denom_vwet > 0) sum(lag1_vwet & lag0_vwet) / denom_vwet else 0

    # Spell adjustments (keep your current approach)
    ds <- dry.spell.change[m]
    ws <- wet.spell.change[m]
    if (!is.finite(ds) || ds <= 0) ds <- 1
    if (!is.finite(ws) || ws <= 0) ws <- 1

    p01 <- p01 / ds
    p02 <- p02 / ds
    p00 <- 1 - p01 - p02
    if (p00 < 0) { p00 <- 0; s <- p01 + p02; if (s > 0) { p01 <- p01 / s; p02 <- p02 / s } }

    p10 <- p10 / ws
    p11 <- 1 - p10 - p12
    if (p11 < 0) { p11 <- 0; s <- p10 + p12; if (s > 0) { p10 <- p10 / s; p12 <- p12 / s } }

    # Assign to all simulation days in that month/year
    p00_final[r] <- p00; p01_final[r] <- p01; p02_final[r] <- p02
    p10_final[r] <- p10; p11_final[r] <- p11; p12_final[r] <- p12
    p20_final[r] <- p20; p21_final[r] <- p21; p22_final[r] <- p22
  }

  list(
    p00_final = p00_final, p01_final = p01_final, p02_final = p02_final,
    p10_final = p10_final, p11_final = p11_final, p12_final = p12_final,
    p20_final = p20_final, p21_final = p21_final, p22_final = p22_final
  )
}

