#' Estimate Monthly Occurrence-State Markov Transition Probabilities
#'
#' @description
#' Estimates month-specific transition probabilities for a three-state
#' precipitation occurrence Markov chain (dry, wet, very wet) from observed
#' daily precipitation data, and maps these probabilities onto a simulated
#' time axis.
#'
#' Transition probabilities are inferred from observed lag-1 to lag-0
#' precipitation pairs using month-dependent wet and extreme thresholds.
#' A Dirichlet-type smoothing is applied to avoid zero-probability transitions,
#' with the effective smoothing strength decreasing automatically with the
#' number of observed transitions per month.
#'
#' The function supports both calendar-year and water-year simulations.
#' The regime is inferred from month_list: if month_list[1] == 1,
#' calendar-year logic is used; otherwise, water-year logic is assumed.
#'
#' Optional spell-persistence modifiers are applied to dry and wet states to
#' influence transition persistence. After these adjustments, transition
#' probabilities are explicitly normalized to ensure stochastic validity
#' (non-negativity and rows summing to one).
#'
#' The estimated monthly transition probabilities are assigned to all simulated
#' days that fall within the target simulation year and corresponding calendar
#' month.
#'
#' @param PRCP_LAG0 Numeric vector. Observed daily precipitation at lag 0
#'   (current day) used to estimate state transitions.
#' @param PRCP_LAG1 Numeric vector. Observed daily precipitation at lag 1
#'   (previous day) used to estimate state transitions.
#' @param MONTH_LAG0 Integer vector (1-12). Calendar month corresponding to
#'   PRCP_LAG0.
#' @param MONTH_LAG1 Integer vector (1-12). Calendar month corresponding to
#'   PRCP_LAG1.
#' @param YEAR_LAG0 Optional integer vector. Calendar year or water year
#'   corresponding to PRCP_LAG0. If provided together with YEAR_LAG1,
#'   transitions are restricted to same-year pairs.
#' @param YEAR_LAG1 Optional integer vector. Calendar year or water year
#'   corresponding to PRCP_LAG1.
#' @param wet_threshold Numeric vector of length 12. Monthly precipitation
#'   thresholds separating dry and wet states, aligned to month_list.
#' @param extreme_threshold Numeric vector of length 12. Monthly precipitation
#'   thresholds separating wet and very wet states, aligned to month_list.
#' @param month_list Integer vector of length 12 defining the simulated ordering
#'   of months (for example 1:12 for calendar years or c(10:12, 1:9) for
#'   October-start water years).
#' @param MONTH_SIM Integer vector. Calendar month for each simulated day.
#' @param WATER_YEAR_SIM Integer vector. Simulation year (calendar or water year,
#'   depending on regime) for each simulated day.
#' @param y Integer. Index of the current simulated year (1-based).
#' @param START_YEAR_SIM Integer. First year of the simulation, interpreted as a
#'   calendar year or water year depending on the simulation regime.
#' @param dry.spell.change Numeric vector of length 12. Monthly multiplicative
#'   factors controlling dry-state persistence.
#' @param wet.spell.change Numeric vector of length 12. Monthly multiplicative
#'   factors controlling wet-state persistence.
#' @param SIM_LENGTH Integer. Total number of simulated days. Determines the
#'   length of the returned probability vectors.
#' @param alpha Numeric. Base smoothing strength for transition-probability
#'   estimation. The effective smoothing applied in each month is
#'   alpha / sqrt(N_m), where N_m is the number of observed transitions
#'   in that month.
#'
#' @return
#' A named list of nine numeric vectors, each of length SIM_LENGTH,
#' representing time-varying transition probabilities:
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
#' Precipitation occurrence states are defined as:
#' \describe{
#'   \item{Dry}{Precipitation less than or equal to the monthly wet threshold}
#'   \item{Wet}{Precipitation greater than the wet threshold and less than or
#'              equal to the extreme threshold}
#'   \item{Very wet}{Precipitation greater than the extreme threshold}
#' }
#'
#' All transition-probability rows are explicitly normalized after smoothing and
#' spell adjustments to ensure valid Markov-chain behavior.
#'
#' @export

monthly_markov_probs <- function(
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
    SIM_LENGTH,
    alpha = 1.0
) {

  if (!is.finite(alpha) || alpha < 0) {
    stop("alpha must be a non-negative finite number")
  }

  # Helper function
  normalize_probs <- function(p, fallback = NULL) {
    p[!is.finite(p) | p < 0] <- 0
    s <- sum(p)
    if (s > 0) {
      p / s
    } else {
      if (is.null(fallback)) {
        rep(1 / length(p), length(p))
      } else {
        fallback
      }
    }
  }

  water.year <- (month_list[1] != 1)
  target_year <- if (water.year) (START_YEAR_SIM + y) else (START_YEAR_SIM + y - 1)

  # Check if spell adjustments should be applied
  use_spell_adjust <- any(abs(dry.spell.change - 1) > 1e-10) ||
    any(abs(wet.spell.change - 1) > 1e-10)

  # Strong lag guard
  if (!is.null(YEAR_LAG0) && !is.null(YEAR_LAG1)) {

    keep <- YEAR_LAG0 == YEAR_LAG1
    PRCP_LAG0  <- PRCP_LAG0[keep]
    PRCP_LAG1  <- PRCP_LAG1[keep]
    MONTH_LAG0 <- MONTH_LAG0[keep]
    MONTH_LAG1 <- MONTH_LAG1[keep]

  } else if (!water.year) {

    keep <- MONTH_LAG0 >= MONTH_LAG1
    PRCP_LAG0  <- PRCP_LAG0[keep]
    PRCP_LAG1  <- PRCP_LAG1[keep]
    MONTH_LAG0 <- MONTH_LAG0[keep]
    MONTH_LAG1 <- MONTH_LAG1[keep]
  }

  # Vectorized state classification

  # Map each observation to month index
  idx_lag1 <- match(MONTH_LAG1, month_list)
  idx_lag0 <- match(MONTH_LAG0, month_list)

  # Get thresholds for each observation
  wet_thr_lag1 <- wet_threshold[idx_lag1]
  ext_thr_lag1 <- extreme_threshold[idx_lag1]
  wet_thr_lag0 <- wet_threshold[idx_lag0]
  ext_thr_lag0 <- extreme_threshold[idx_lag0]

  # Handle non-finite thresholds (fallback to global)
  if (any(!is.finite(wet_thr_lag1)) || any(!is.finite(ext_thr_lag1))) {
    global_wet <- quantile(PRCP_LAG1, 0.2, na.rm = TRUE)
    pos <- PRCP_LAG1[PRCP_LAG1 > 0]
    global_ext <- quantile(if(length(pos)) pos else PRCP_LAG1, 0.8, na.rm = TRUE)

    wet_thr_lag1[!is.finite(wet_thr_lag1)] <- global_wet
    ext_thr_lag1[!is.finite(ext_thr_lag1)] <- global_ext
    wet_thr_lag0[!is.finite(wet_thr_lag0)] <- global_wet
    ext_thr_lag0[!is.finite(ext_thr_lag0)] <- global_ext
  }

  # Classify states (vectorized - single pass)
  state_lag1 <- ifelse(PRCP_LAG1 <= wet_thr_lag1, 0L,
                       ifelse(PRCP_LAG1 <= ext_thr_lag1, 1L, 2L))

  state_lag0 <- ifelse(PRCP_LAG0 <= wet_thr_lag0, 0L,
                       ifelse(PRCP_LAG0 <= ext_thr_lag0, 1L, 2L))

  ## -------------------------------------------------
  ## Output containers
  ## -------------------------------------------------
  p00_final <- p01_final <- p02_final <- rep(NA_real_, SIM_LENGTH)
  p10_final <- p11_final <- p12_final <- rep(NA_real_, SIM_LENGTH)
  p20_final <- p21_final <- p22_final <- rep(NA_real_, SIM_LENGTH)

  ## -------------------------------------------------
  ## Month loop
  ## -------------------------------------------------
  for (m in seq_along(month_list)) {

    mm <- month_list[m]
    r <- which(MONTH_SIM == mm & WATER_YEAR_SIM == target_year)
    if (!length(r)) next

    x <- which(MONTH_LAG1 == mm)
    N_m <- length(x)
    alpha_m <- if (N_m > 0L) alpha / sqrt(N_m) else alpha

    if (!length(x)) {
      p00_final[r] <- 1; p01_final[r] <- 0; p02_final[r] <- 0
      p10_final[r] <- 1; p11_final[r] <- 0; p12_final[r] <- 0
      p20_final[r] <- 1; p21_final[r] <- 0; p22_final[r] <- 0
      next
    }

    s1 <- state_lag1[x]
    s0 <- state_lag0[x]

    ## ---- Raw transition counts ----
    n00 <- sum(s1 == 0 & s0 == 0)
    n01 <- sum(s1 == 0 & s0 == 1)
    n02 <- sum(s1 == 0 & s0 == 2)

    n10 <- sum(s1 == 1 & s0 == 0)
    n11 <- sum(s1 == 1 & s0 == 1)
    n12 <- sum(s1 == 1 & s0 == 2)

    n20 <- sum(s1 == 2 & s0 == 0)
    n21 <- sum(s1 == 2 & s0 == 1)
    n22 <- sum(s1 == 2 & s0 == 2)

    ## ---- Dirichlet smoothing (applied BEFORE spell adjustment) ----
    p_dry <- c(n00, n01, n02) + alpha_m
    p_wet <- c(n10, n11, n12) + alpha_m
    p_vwt <- c(n20, n21, n22) + alpha_m

    p_dry <- p_dry / sum(p_dry)
    p_wet <- p_wet / sum(p_wet)
    p_vwt <- p_vwt / sum(p_vwt)

    ## -------------------------------------------------
    ## SPELL ADJUSTMENTS (conditional on user input)
    ## -------------------------------------------------
    if (use_spell_adjust) {
      ds <- dry.spell.change[m]
      ws <- wet.spell.change[m]

      # Validate factors for this month
      if (!is.finite(ds) || ds <= 0) ds <- 1
      if (!is.finite(ws) || ws <= 0) ws <- 1

      # Apply dry spell adjustment only if factor differs from 1
      if (abs(ds - 1) > 1e-10) {
        # Reduce transitions away from dry state
        p_dry <- normalize_probs(
          c(p_dry[1], p_dry[2] / ds, p_dry[3] / ds),
          fallback = c(1, 0, 0)
        )
      }

      # Apply wet spell adjustment only if factor differs from 1
      if (abs(ws - 1) > 1e-10) {
        # Reduce transition to dry state from wet
        p_wet <- normalize_probs(
          c(p_wet[1] / ws, p_wet[2], p_wet[3]),
          fallback = c(0, 1, 0)
        )
      }

      # Very wet row: no adjustment applied (already normalized)
    }

    ## ---- Assign probabilities to output vectors ----
    # CRITICAL: This must be OUTSIDE the spell adjustment if block
    # so probabilities are assigned whether or not adjustments are applied
    p00_final[r] <- p_dry[1]; p01_final[r] <- p_dry[2]; p02_final[r] <- p_dry[3]
    p10_final[r] <- p_wet[1]; p11_final[r] <- p_wet[2]; p12_final[r] <- p_wet[3]
    p20_final[r] <- p_vwt[1]; p21_final[r] <- p_vwt[2]; p22_final[r] <- p_vwt[3]
  }

  list(
    p00_final = p00_final, p01_final = p01_final, p02_final = p02_final,
    p10_final = p10_final, p11_final = p11_final, p12_final = p12_final,
    p20_final = p20_final, p21_final = p21_final, p22_final = p22_final
  )
}

