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
#' The regime is inferred from month.list: if month.list[1] == 1,
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
#' @param precip.lag0 Numeric vector. Observed daily precipitation at lag 0
#'   (current day) used to estimate state transitions.
#' @param precip.lag1 Numeric vector. Observed daily precipitation at lag 1
#'   (previous day) used to estimate state transitions.
#' @param month.lag0 Integer vector (1-12). Calendar month corresponding to
#'   precip.lag0.
#' @param month.lag1 Integer vector (1-12). Calendar month corresponding to
#'   precip.lag1.
#' @param year.lag0 Optional integer vector. Calendar year or water year
#'   corresponding to precip.lag0. If provided together with year.lag1,
#'   transitions are restricted to same-year pairs.
#' @param year.lag1 Optional integer vector. Calendar year or water year
#'   corresponding to precip.lag1.
#' @param wet.threshold Numeric vector of length 12. Monthly precipitation
#'   thresholds separating dry and wet states, aligned to month.list.
#' @param extreme.threshold Numeric vector of length 12. Monthly precipitation
#'   thresholds separating wet and very wet states, aligned to month.list.
#' @param month.list Integer vector of length 12 defining the simulated ordering
#'   of months (for example 1:12 for calendar years or c(10:12, 1:9) for
#'   October-start water years).
#' @param sim.months Integer vector. Calendar month for each simulated day.
#' @param sim.water.years Integer vector. Simulation year (calendar or water year,
#'   depending on regime) for each simulated day.
#' @param year.idx Integer. Index of the current simulated year (1-based).
#' @param sim.start.year Integer. First year of the simulation, interpreted as a
#'   calendar year or water year depending on the simulation regime.
#' @param dry.spell.change Numeric vector of length 12. Monthly multiplicative
#'   factors controlling dry-state persistence.
#' @param wet.spell.change Numeric vector of length 12. Monthly multiplicative
#'   factors controlling wet-state persistence.
#' @param sim.length Integer. Total number of simulated days. Determines the
#'   length of the returned probability vectors.
#' @param alpha Numeric. Base smoothing strength for transition-probability
#'   estimation. The effective smoothing applied in each month is
#'   alpha / sqrt(N_m), where N_m is the number of observed transitions
#'   in that month.
#'
#' @return
#' A named list of nine numeric vectors, each of length sim.length,
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
    precip.lag0,
    precip.lag1,
    month.lag0,
    month.lag1,
    year.lag0 = NULL,
    year.lag1 = NULL,
    wet.threshold,
    extreme.threshold,
    month.list,
    sim.months,
    sim.water.years,
    year.idx,
    sim.start.year,
    dry.spell.change,
    wet.spell.change,
    sim.length,
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

  water.year <- (month.list[1] != 1)
  target_year <- if (water.year) (sim.start.year + year.idx) else (sim.start.year + year.idx - 1)

  # Check if spell adjustments should be applied
  use_spell_adjust <- any(abs(dry.spell.change - 1) > 1e-10) ||
    any(abs(wet.spell.change - 1) > 1e-10)

  # Strong lag guard
  if (!is.null(year.lag0) && !is.null(year.lag1)) {

    keep <- year.lag0 == year.lag1
    precip.lag0  <- precip.lag0[keep]
    precip.lag1  <- precip.lag1[keep]
    month.lag0 <- month.lag0[keep]
    month.lag1 <- month.lag1[keep]

  } else if (!water.year) {

    keep <- month.lag0 >= month.lag1
    precip.lag0  <- precip.lag0[keep]
    precip.lag1  <- precip.lag1[keep]
    month.lag0 <- month.lag0[keep]
    month.lag1 <- month.lag1[keep]
  }

  # Vectorized state classification

  # Map each observation to month index
  idx_lag1 <- match(month.lag1, month.list)
  idx_lag0 <- match(month.lag0, month.list)

  # Get thresholds for each observation
  wet_thr_lag1 <- wet.threshold[idx_lag1]
  ext_thr_lag1 <- extreme.threshold[idx_lag1]
  wet_thr_lag0 <- wet.threshold[idx_lag0]
  ext_thr_lag0 <- extreme.threshold[idx_lag0]

  # Handle non-finite thresholds (fallback to global)
  if (any(!is.finite(wet_thr_lag1)) || any(!is.finite(ext_thr_lag1))) {
    global_wet <- quantile(precip.lag1, 0.2, na.rm = TRUE)
    pos <- precip.lag1[precip.lag1 > 0]
    global_ext <- quantile(if(length(pos)) pos else precip.lag1, 0.8, na.rm = TRUE)

    wet_thr_lag1[!is.finite(wet_thr_lag1)] <- global_wet
    ext_thr_lag1[!is.finite(ext_thr_lag1)] <- global_ext
    wet_thr_lag0[!is.finite(wet_thr_lag0)] <- global_wet
    ext_thr_lag0[!is.finite(ext_thr_lag0)] <- global_ext
  }

  # Classify states (vectorized - single pass)
  state_lag1 <- ifelse(precip.lag1 <= wet_thr_lag1, 0L,
                       ifelse(precip.lag1 <= ext_thr_lag1, 1L, 2L))

  state_lag0 <- ifelse(precip.lag0 <= wet_thr_lag0, 0L,
                       ifelse(precip.lag0 <= ext_thr_lag0, 1L, 2L))

  ## -------------------------------------------------
  ## Output containers
  ## -------------------------------------------------
  p00_final <- p01_final <- p02_final <- rep(NA_real_, sim.length)
  p10_final <- p11_final <- p12_final <- rep(NA_real_, sim.length)
  p20_final <- p21_final <- p22_final <- rep(NA_real_, sim.length)

  ## -------------------------------------------------
  ## Month loop
  ## -------------------------------------------------
  for (m in seq_along(month.list)) {

    mm <- month.list[m]
    r <- which(sim.months == mm & sim.water.years == target_year)
    if (!length(r)) next

    x <- which(month.lag1 == mm)
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




#' Determine Next State in a First-Order Markov Chain
#'
#' Given the previous state and a random number `rn`, this function returns the next
#' state in a 3-state first-order Markov chain. The transition probabilities are
#' state- and index-specific.
#'
#' States are typically defined as:
#' - 0: Dry
#' - 1: Wet
#' - 2: Extreme
#'
#' The transition is determined using cumulative probabilities for the transitions
#' to state 0 and 1. Any remaining probability is assigned to state 2.
#'
#' @param prev_state Integer. The current (previous) state (0, 1, or 2).
#' @param rn Numeric. A random number in [0, 1] used to sample the next state.
#' @param idx Integer. Index for selecting transition probabilities (e.g., time step or spatial location).
#' @param p00 Numeric vector. Probability of transition from state 0 to 0.
#' @param p01 Numeric vector. Probability of transition from state 0 to 1.
#' @param p10 Numeric vector. Probability of transition from state 1 to 0.
#' @param p11 Numeric vector. Probability of transition from state 1 to 1.
#' @param p20 Numeric vector. Probability of transition from state 2 to 0.
#' @param p21 Numeric vector. Probability of transition from state 2 to 1.
#'
#' @return Integer. The next state (0, 1, or 2).
#'
#' @examples
#' # Suppose we're at time step 5 with current state = 1
#' set.seed(123)
#' rn <- runif(1)
#' markov_next_state(
#'   prev_state = 1,
#'   rn = rn,
#'   idx = 5,
#'   p00 = rep(0.7, 10),
#'   p01 = rep(0.2, 10),
#'   p10 = rep(0.3, 10),
#'   p11 = rep(0.4, 10),
#'   p20 = rep(0.1, 10),
#'   p21 = rep(0.3, 10)
#' )
#'
#' @export
markov_next_state <- function(prev_state, rn, idx, p00, p01, p10, p11, p20, p21) {

  # Clamp random number
  if (!is.finite(rn)) rn <- 0
  if (rn < 0) rn <- 0
  if (rn > 1) rn <- 1

  # Validate prev_state
  if (!(prev_state %in% c(0L, 1L, 2L))) {
    prev_state <- 0L
  }

  # Clamp idx to valid range
  n <- length(p00)
  if (idx < 1L) idx <- 1L
  if (idx > n)  idx <- n

  # Select row probabilities
  if (prev_state == 0L) {
    a <- p00[idx]
    b <- p01[idx]
  } else if (prev_state == 1L) {
    a <- p10[idx]
    b <- p11[idx]
  } else {
    a <- p20[idx]
    b <- p21[idx]
  }

  # NA -> 0
  if (!is.finite(a)) a <- 0
  if (!is.finite(b)) b <- 0

  # Clamp negatives
  if (a < 0) a <- 0
  if (b < 0) b <- 0

  # Enforce stochastic validity
  s <- a + b
  if (s > 1) {
    if (s > 1.01) {  # Tolerance for rounding
      warning(sprintf("Invalid transition probabilities: sum = %.3f at idx %d", s, idx))
    }
    a <- a / s
    b <- b / s
    s <- 1  # Force exact normalization
  }

  # Draw next state
  if (rn < a) {
    0L
  } else if (rn < (a + b)) {
    1L
  } else {
    2L
  }
}

#' Get Indices for a Specific Occurrence State Transition
#'
#' Returns the indices in a precipitation time series (prcp) where a specified
#' state transition occurs, based on wet/dry/extreme day thresholds.
#'
#' This optimized version pre-classifies all states once, avoiding repeated
#' subsetting and threshold comparisons. Typically 5-10x faster than the
#' if-else chain approach.
#'
#' The occurrence states are defined as:
#' - 0: Dry day (precipitation <= wet.thr)
#' - 1: Wet day (precipitation > wet.thr and <= extreme.thr)
#' - 2: Extreme wet day (precipitation > extreme.thr)
#'
#' @param from.state Integer. The current occurrence state (0 = dry, 1 = wet, 2 = extreme).
#' @param to.state Integer. The next occurrence state (0 = dry, 1 = wet, 2 = extreme).
#' @param prcp Numeric vector of precipitation values.
#' @param candidate.idx Integer vector of indices representing the current day in the time series.
#' @param wet.thr Numeric. Threshold separating dry and wet days.
#' @param extreme.thr Numeric. Threshold above which days are considered extreme.
#'
#' @return Integer vector of indices where the given transition from from.state to to.state occurs.
#'
#' @details
#' This function pre-extracts and classifies precipitation states once, then
#' performs a simple integer comparison to find matching transitions. This is
#' much more efficient than the if-else chain approach which would subset and
#' compare thresholds multiple times.
#'
#' @examples
#' prcp <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
#' candidate.idx <- 1:(length(prcp) - 1)
#' wet.thr <- 1
#' extreme.thr <- 20
#'
#' # Find transitions from dry to wet
#' get_state_indices(0, 1, prcp, candidate.idx, wet.thr, extreme.thr)
#' # Returns: 6 (index where transition occurs)
#'
#' # Find transitions from wet to extreme
#' get_state_indices(1, 2, prcp, candidate.idx, wet.thr, extreme.thr)
#' # Returns: 3 7 (indices where transitions occur)
#'
#' # Find transitions from extreme to dry
#' get_state_indices(2, 0, prcp, candidate.idx, wet.thr, extreme.thr)
#' # Returns: 5 9 (indices where transitions occur)
#'
#' @export
get_state_indices <- function(from.state, to.state, prcp, candidate.idx, wet.thr, extreme.thr) {

  #Extract precipitation values
  prcp_cur <- prcp[candidate.idx]
  prcp_next <- prcp[candidate.idx + 1]

  state_cur <- ifelse(prcp_cur <= wet.thr, 0L,
                      ifelse(prcp_cur <= extreme.thr, 1L, 2L))

  state_next <- ifelse(prcp_next <= wet.thr, 0L,
                       ifelse(prcp_next <= extreme.thr, 1L, 2L))

  # Returns indices in candidate.idx where the transition occurs
  which(state_cur == from.state & state_next == to.state)
}


#' Safely Get or Sample an Index from a Vector
#'
#' Returns a valid index from `precip.tomorrow` based on the provided `result`.
#' If `result` is missing, out of bounds, or invalid, a random index is sampled instead.
#' If `precip.tomorrow` is empty, `NA_integer_` is returned.
#'
#' This function is useful in stochastic weather generation workflows where fallback
#' behavior is needed when a computed index is not valid.
#'
#' @param result Integer. Proposed index for retrieving an element from `precip.tomorrow`.
#' @param precip.tomorrow Numeric vector. A set of candidate precipitation values (e.g., for the next day).
#'
#' @return Integer. A valid index from `precip.tomorrow`, or `NA_integer_` if `precip.tomorrow` is empty.
#'
#' @examples
#' precip.tomorrow <- c(0, 5, 10, 20)
#'
#' # Valid index
#' get_result_index(2, precip.tomorrow) # returns 2
#'
#' # Invalid index: falls back to sampling
#' set.seed(1)
#' get_result_index(10, precip.tomorrow) # randomly returns a valid index
#'
#' # NA index: falls back to sampling
#' get_result_index(NA, precip.tomorrow)
#'
#' # Empty vector: returns NA
#' get_result_index(1, numeric(0))
#'
#' @export
get_result_index <- function(result, precip.tomorrow) {
  if (is.na(result) || length(result) == 0 || result < 1 || result > length(precip.tomorrow)) {
    if (length(precip.tomorrow) > 0) {
      sample(seq_along(precip.tomorrow), 1)
    } else {
      NA_integer_
    }
  } else {
    result
  }
}
