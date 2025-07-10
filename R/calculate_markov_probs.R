#' Calculate Monthly Markov Transition Probabilities (Vectorized)
#'
#' @description
#' Computes the monthly transition probabilities for a three-state Markov chain (dry, wet, very wet),
#' using daily lagged precipitation data and monthly thresholds. This function is fully vectorized for efficiency
#' and is designed for use within daily stochastic weather generators.
#'
#' @param PRCP_LAG0 Numeric vector, lag 0 (current day) precipitation for all days.
#' @param PRCP_LAG1 Numeric vector, lag 1 (previous day) precipitation for all days.
#' @param MONTH_LAG0 Integer vector, month (1-12) for lag 0 (current day) for all days.
#' @param MONTH_LAG1 Integer vector, month (1-12) for lag 1 (previous day) for all days.
#' @param wet_threshold Numeric vector of length 12, monthly thresholds separating dry/wet.
#' @param extreme_threshold Numeric vector of length 12, monthly thresholds for "very wet" state.
#' @param month_list Integer vector of length 12, mapping simulation months to calendar months.
#' @param MONTH_SIM Integer vector, simulation months for each simulated day.
#' @param WATER_YEAR_SIM Integer vector, simulation water year for each simulated day.
#' @param y Integer, current simulated year index.
#' @param START_YEAR_SIM Integer, simulation's starting year.
#' @param dry.spell.change Numeric vector of length 12, adjustment factors for dry spell persistence.
#' @param wet.spell.change Numeric vector of length 12, adjustment factors for wet spell persistence.
#' @param SIM_LENGTH Integer, length of the simulation time series.
#'
#' @return A list with nine elements, each a numeric vector of length `SIM_LENGTH`:
#'   \item{p00_final}{P(dry???dry)}
#'   \item{p01_final}{P(dry???wet)}
#'   \item{p02_final}{P(dry???very wet)}
#'   \item{p10_final}{P(wet???dry)}
#'   \item{p11_final}{P(wet???wet)}
#'   \item{p12_final}{P(wet???very wet)}
#'   \item{p20_final}{P(very wet???dry)}
#'   \item{p21_final}{P(very wet???wet)}
#'   \item{p22_final}{P(very wet???very wet)}
#'
#' @details
#' The function computes monthly state transition probabilities for a 3-state weather occurrence Markov chain,
#' with states defined by precipitation thresholds (dry, wet, very wet). Probabilities are adjusted
#' for spell persistence via user-supplied adjustment vectors.
#'
#' @export
calculateMarkovProbs <- function(PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
                                 wet_threshold, extreme_threshold,
                                 month_list, MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
                                 dry.spell.change, wet.spell.change, SIM_LENGTH) {
  # Preallocate storage
  p00_final <- numeric(SIM_LENGTH)
  p01_final <- numeric(SIM_LENGTH)
  p02_final <- numeric(SIM_LENGTH)
  p10_final <- numeric(SIM_LENGTH)
  p11_final <- numeric(SIM_LENGTH)
  p12_final <- numeric(SIM_LENGTH)
  p20_final <- numeric(SIM_LENGTH)
  p21_final <- numeric(SIM_LENGTH)
  p22_final <- numeric(SIM_LENGTH)

  for (m in 1:12) {
    x <- which(MONTH_LAG1 == month_list[m])
    r <- which(MONTH_SIM == month_list[m] & WATER_YEAR_SIM == (y + START_YEAR_SIM))
    if (length(x) == 0 || length(r) == 0) next

    # Vectorized state logic (no nested for-loops)
    lag1_wet <- PRCP_LAG1[x] > wet_threshold[m] & PRCP_LAG1[x] <= extreme_threshold[m]
    lag1_dry <- PRCP_LAG1[x] <= wet_threshold[m]
    lag1_vwet <- PRCP_LAG1[x] > extreme_threshold[m]
    lag0_dry <- PRCP_LAG0[x] <= wet_threshold[m]
    lag0_wet <- PRCP_LAG0[x] > wet_threshold[m] & PRCP_LAG0[x] <= extreme_threshold[m]
    lag0_vwet <- PRCP_LAG0[x] > extreme_threshold[m]

    denom_dry <- sum(lag1_dry)
    denom_wet <- sum(lag1_wet)
    denom_vwet <- sum(lag1_vwet)

    # Calculate each transition probability
    p00_final[r] <- if (denom_dry > 0) sum(lag1_dry & lag0_dry) / denom_dry else 0
    p01_final[r] <- if (denom_dry > 0) sum(lag1_dry & lag0_wet) / denom_dry else 0
    p02_final[r] <- if (denom_dry > 0) sum(lag1_dry & lag0_vwet) / denom_dry else 0
    p10_final[r] <- if (denom_wet > 0) sum(lag1_wet & lag0_dry) / denom_wet else 0
    p11_final[r] <- if (denom_wet > 0) sum(lag1_wet & lag0_wet) / denom_wet else 0
    p12_final[r] <- if (denom_wet > 0) sum(lag1_wet & lag0_vwet) / denom_wet else 0
    p20_final[r] <- if (denom_vwet > 0) sum(lag1_vwet & lag0_dry) / denom_vwet else 0
    p21_final[r] <- if (denom_vwet > 0) sum(lag1_vwet & lag0_wet) / denom_vwet else 0
    p22_final[r] <- if (denom_vwet > 0) sum(lag1_vwet & lag0_vwet) / denom_vwet else 0

    # Adjust for dry/wet spell persistence
    p01_final[r] <- p01_final[r] / dry.spell.change[m]
    p02_final[r] <- p02_final[r] / dry.spell.change[m]
    p00_final[r] <- 1 - p01_final[r] - p02_final[r]

    p10_final[r] <- p10_final[r] / wet.spell.change[m]
    p11_final[r] <- 1 - p10_final[r] - p12_final[r]
  }

  list(
    p00_final = p00_final,
    p01_final = p01_final,
    p02_final = p02_final,
    p10_final = p10_final,
    p11_final = p11_final,
    p12_final = p12_final,
    p20_final = p20_final,
    p21_final = p21_final,
    p22_final = p22_final
  )
}
