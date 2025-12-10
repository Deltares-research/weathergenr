#' @title Resample Daily Dates for Stochastic Weather Generation
#'
#' @description
#' Generates a simulated daily weather date sequence using a Markov Chain model and K-Nearest Neighbors (KNN) resampling,
#' conditioned on annual precipitation and climate statistics. This function is intended for stochastic weather generators,
#' producing daily dates that preserve the statistical characteristics and persistence of observed wet/dry/extreme spells.
#'
#' @param PRCP_FINAL_ANNUAL_SIM Numeric vector. Simulated annual precipitation totals, one per simulation year.
#' @param ANNUAL_PRCP Numeric vector. Observed historical annual precipitation totals, one per observed year.
#' @param PRCP Numeric vector. Observed daily precipitation series (length must match `dates.d`).
#' @param TEMP Numeric vector. Observed daily temperature series (e.g., mean temperature; length must match `dates.d`).
#' @param START_YEAR_SIM Integer. First year of the simulation period.
#' @param k1 Integer. Simulation trace index (used as a seed modifier for reproducibility).
#' @param ymax Integer. Number of years to simulate (simulation period length).
#' @param dates.d Data frame. Observed period with columns: `date`, `month`, `day`, `wyear` (water year).
#' @param sim.dates.d Data frame. Simulation period with columns: `month`, `day`, `wyear`.
#' @param YEAR_D Integer vector. Water year for each observed day (length must match `dates.d`).
#' @param month.start Integer (1-12). First month of the water year (e.g., 1 for January).
#' @param knn.annual.sample.num Integer. Number of nearest neighbors to sample for annual KNN matching (default: 50).
#' @param wet.quantile Numeric (0-1). Quantile threshold for defining wet days (default: 0.2).
#' @param extreme.quantile Numeric (0-1). Quantile threshold for defining extreme wet days (default: 0.8).
#' @param dry.spell.change Numeric vector (length 12). Monthly adjustment factors for dry spell lengths (default: rep(1, 12)).
#' @param wet.spell.change Numeric vector (length 12). Monthly adjustment factors for wet spell lengths (default: rep(1, 12)).
#' @param seed Optional integer. Seed for random number generation to ensure reproducibility.
#'
#' @details
#' This function resamples daily weather dates to create a synthetic weather sequence that matches both annual
#' and daily precipitation/temperature statistics of the observed record. The process involves:
#'
#' 1. **Annual selection:** For each simulation year, similar observed years are identified using KNN based on annual precipitation.
#' 2. **Markov Chain simulation:** Daily precipitation states (dry, wet, extreme) are generated using a Markov Chain model, with transition probabilities estimated from the observed years and adjusted for each month.
#' 3. **Daily KNN resampling:** For each simulated day, a corresponding observed day is selected via KNN based on standardized anomalies of precipitation and temperature, to match the simulated state and maintain realistic temporal structure.
#' 4. **Spell adjustments:** `dry.spell.change` and `wet.spell.change` allow monthly control over the lengths of dry and wet spells in the synthetic sequence.
#'
#' The result is a vector of `Date` objects for the synthetic daily weather sequence, preserving the statistical
#' properties of both annual and sub-annual (spell and state) variability.
#'
#' @return
#' A vector of `Date` objects representing the resampled daily dates for the simulation period.
#'
#' @examples
#' \dontrun{
#' sim_dates <- resampleDates(
#'   PRCP_FINAL_ANNUAL_SIM = runif(3, 900, 1100),
#'   ANNUAL_PRCP = runif(20, 900, 1100),
#'   PRCP = rnorm(7300, 3, 5),
#'   TEMP = rnorm(7300, 18, 5),
#'   START_YEAR_SIM = 2001,
#'   k1 = 1,
#'   ymax = 3,
#'   dates.d = data.frame(
#'     date = seq.Date(as.Date("2000-01-01"), as.Date("2019-12-31"), by = "day"),
#'     month = rep(1:12, each = 365, length.out = 7300),
#'     day = rep(1:31, length.out = 7300),
#'     wyear = rep(2000:2019, each = 365)
#'   ),
#'   sim.dates.d = data.frame(
#'     month = rep(1:12, each = 365, length.out = 1095),
#'     day = rep(1:31, length.out = 1095),
#'     wyear = rep(2001:2003, each = 365)
#'   ),
#'   YEAR_D = rep(2000:2019, each = 365),
#'   month.start = 1,
#'   seed = 123
#' )
#' head(sim_dates)
#' }
#'
#' @export
resampleDates <- function(
    PRCP_FINAL_ANNUAL_SIM,
    ANNUAL_PRCP,
    PRCP,
    TEMP,
    START_YEAR_SIM, #redundant
    k1, #? check redundancy
    ymax, #
    dates.d,
    sim.dates.d,
    YEAR_D, # redundant
    month.start = 1,
    knn.annual.sample.num = 50,
    wet.quantile = 0.2,
    extreme.quantile = 0.8,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    seed = NULL

) {


  # -- RNG seed for reproducibility, local only
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed + k1)
  }


  # -- Month order for water years
  month_list <- if (month.start == 1) 1:12 else c(month.start:12, 1:(month.start - 1))

  # -- Set observed date vectors
  dates_obs <- dates.d$date
  months_obs <- dates.d$month
  wyears_obs <- dates.d$wyear
  monthday_obs <- dates.d[, c("month", "day")]
  water_years_obs <- dates.d$wyear[dates.d$month == month.start & dates.d$day == 1]

  # -- Pre-compute observed weather statistics
  mean_mon_TEMP <- tapply(TEMP, months_obs, mean)
  mean_mon_PRCP <- tapply(PRCP, months_obs, mean)
  sd_mon_TEMP <- tapply(TEMP, months_obs, sd)
  sd_mon_PRCP <- tapply(PRCP, months_obs, sd)


  # -- Set simulated date vectors
  sim_length <- nrow(sim.dates.d)
  sim_date <- as.Date(rep(NA, sim_length))
  months_sim <- sim.dates.d$month
  days_sim <- sim.dates.d$day
  wyears_sim <- sim.dates.d$wyear


  # --- Initialize first simulated day's value sensibly ---
  first_month <- months_sim[1]
  first_day <- days_sim[1]

  # Use all obs with the same month/day as the first simulated day
  first_obs_idx <- which(months_obs == first_month & monthday_obs$day == first_day)
  if (length(first_obs_idx) == 0) first_obs_idx <- 1 # fallback: just use first obs if no match

  # -- Set simulation vector and initial states
  sim_prcp <- numeric(sim_length)
  sim_temp <- numeric(sim_length)
  sim_occurrence <- integer(sim_length)

  sim_prcp[1] <- PRCP[sample(first_obs_idx, 1)]
  sim_temp[1] <- TEMP[sample(first_obs_idx, 1)]
  sim_date[1] <- as.Date(dates_obs[sample(first_obs_idx, 1)])
  sim_occurrence[1] <- 0 # dry state by default, or you can sample from observed state if desired

  # -- Random draws upfront
  rn_all <- stats::runif(sim_length)

  # -- k for annual KNN
  k_annual <- max(round(sqrt(length(ANNUAL_PRCP))), round(length(ANNUAL_PRCP) * 0.5))

  # -- Precompute daily lookup for observed days by month-day (for fast subsetting)
  monthday_key <- paste(monthday_obs$month, monthday_obs$day, sep = ".")
  lookup_day_idx <- split(seq_along(monthday_key), monthday_key)

  sim_idx <- 1

  # ==== Main Simulation Loop (Years) ====
  for (y in seq_len(ymax)) {

    # --- Annual KNN: find years in obs similar to this sim year ---
    year_sample_idx <- knn_sample(
      candidates = ANNUAL_PRCP,
      target = PRCP_FINAL_ANNUAL_SIM[y],
      k = k_annual,
      n = knn.annual.sample.num,
      prob = TRUE,
      weights = NULL,
      seed = seed + k1 * y
    )
    cur_years <- water_years_obs[year_sample_idx]
    obs_idx <- which(wyears_obs %in% cur_years)

    # Subset all vars for current selected years
    prcp_y <- PRCP[obs_idx]
    temp_y <- TEMP[obs_idx]
    date_y <- dates_obs[obs_idx]
    month_y <- months_obs[obs_idx]
    monthday_y <- monthday_obs[obs_idx, ]

    # --- Thresholds per month for current selection ---
    wet_thresh <- rep(stats::quantile(prcp_y, wet.quantile, names = FALSE), 12)
    extreme_thresh <- sapply(seq_along(month_list), function(m) {
      vals <- prcp_y[month_y == month_list[m] & prcp_y > wet_thresh[m]]
      if (length(vals) == 0) NA_real_ else stats::quantile(vals, extreme.quantile, names = FALSE)
    })

    # -- Markov transition probabilities
    probs <- calculateMarkovProbs(
      PRCP_LAG0 = prcp_y[-1],
      PRCP_LAG1 = prcp_y[-length(prcp_y)],
      MONTH_LAG0 = month_y[-1],
      MONTH_LAG1 = month_y[-length(prcp_y)],
      wet_threshold = wet_thresh,
      extreme_threshold = extreme_thresh,
      month_list,
      MONTH_SIM = months_sim,
      WATER_YEAR_SIM = wyears_sim,
      y,
      START_YEAR_SIM,
      dry.spell.change,
      wet.spell.change,
      sim_length
    )

    # Jointly sample initial day's values and state ====
    first_month <- months_sim[sim_idx]
    first_day   <- days_sim[sim_idx]
    mmm <- which(month_list == first_month)
    first_day_indices <- which(month_y == first_month & monthday_y$day == first_day)
    if (length(first_day_indices) == 0) first_day_indices <- 1 # fallback

    prcp_candidates <- prcp_y[first_day_indices]
    temp_candidates <- temp_y[first_day_indices]
    date_candidates <- date_y[first_day_indices]

    # Compute state for each candidate using year-specific thresholds
    candidate_states <- ifelse(
      prcp_candidates == 0, 0,
      ifelse(
        !is.na(extreme_thresh[mmm]) & prcp_candidates > extreme_thresh[mmm], 2,
        ifelse(
          !is.na(wet_thresh[mmm]) & prcp_candidates > wet_thresh[mmm], 1, 0
        )
      )
    )

    # Sample one index
    sample_idx <- sample(seq_along(prcp_candidates), 1)

    sim_prcp[sim_idx]       <- prcp_candidates[sample_idx]
    sim_occurrence[sim_idx] <- candidate_states[sample_idx]
    sim_temp[sim_idx]       <- temp_candidates[sample_idx]
    sim_date[sim_idx]       <- as.Date(date_candidates[sample_idx])


    # -- For each day in the simulation year (assumes 365-day years!)
    for (j in seq_len(365)) {

      sim_idx <- sim_idx + 1
      prev_idx <- sim_idx - 1
      if (sim_idx > sim_length) break

      # --- Markov: state for current day
      sim_occurrence[sim_idx] <- markov_next_state(
        sim_occurrence[prev_idx], rn_all[prev_idx], prev_idx,
        probs$p00_final, probs$p01_final, probs$p10_final,
        probs$p11_final, probs$p20_final, probs$p21_final
      )

      # -- Simulated month/day for current day
      cur_month <- months_sim[prev_idx]
      cur_day <- days_sim[prev_idx]
      m_idx <- which(month_list == cur_month)

      # --- Candidate observed days with same month/day (plus neighbors)
      obs_key <- paste(cur_month, cur_day, sep = ".")
      obs_candidates <- lookup_day_idx[[obs_key]]
      obs_candidates <- obs_candidates[obs_candidates > 0 & obs_candidates <= length(prcp_y)]
      obs_window <- c(-3, -2, -1, 0, 1, 2, 3)
      obs_idx_window <- unique(as.vector(outer(obs_candidates, obs_window, "+")))
      obs_idx_window <- obs_idx_window[obs_idx_window > 0 & obs_idx_window < length(prcp_y)] # fix here

      # Find possible days matching simulated state
      cur_state <- sim_occurrence[prev_idx]
      next_state <- sim_occurrence[sim_idx]
      state_idx <- get_state_indices(
        cur_state, next_state, prcp_y, obs_idx_window,
        wet_thresh[cur_month], extreme_thresh[m_idx]
      )

      # If no matches, expand search window
      if (length(state_idx) == 0) {
        obs_idx_large <- unique(as.vector(outer(obs_candidates, seq(-30, 30), "+")))
        obs_idx_large <- obs_idx_large[obs_idx_large > 0 & obs_idx_large < length(prcp_y)]
        state_idx <- get_state_indices(
          cur_state, next_state, prcp_y, obs_idx_large,
          wet_thresh[cur_month], extreme_thresh[m_idx]
        )
      }

      if (length(state_idx) == 0) {
        # No candidates: fallback to previous day's value
        sim_prcp[sim_idx] <- sim_prcp[prev_idx]
        sim_temp[sim_idx] <- sim_temp[prev_idx]
        sim_date[sim_idx] <- sim_date[prev_idx]
        next
      }

      possible_days <- obs_idx_window[state_idx]
      possible_days <- possible_days[possible_days > 0 & possible_days < length(prcp_y)]

      if (length(possible_days) == 0) {
        sim_prcp[sim_idx] <- sim_prcp[prev_idx]
        sim_temp[sim_idx] <- sim_temp[prev_idx]
        sim_date[sim_idx] <- sim_date[prev_idx]
        next
      }

      # KNN sampling for day selection
      k_day <- max(1, round(sqrt(length(possible_days))))
      prcp_today <- prcp_y[possible_days]
      temp_today <- temp_y[possible_days]

      # Only keep possible_days where possible_days+1 is valid
      valid_idx <- which(possible_days + 1 <= length(prcp_y))
      if (length(valid_idx) == 0) {
        sim_prcp[sim_idx] <- sim_prcp[prev_idx]
        sim_temp[sim_idx] <- sim_temp[prev_idx]
        sim_date[sim_idx] <- sim_date[prev_idx]
        next
      }

      possible_days <- possible_days[valid_idx]
      prcp_today <- prcp_today[valid_idx]
      temp_today <- temp_today[valid_idx]
      prcp_tomorrow <- prcp_y[possible_days + 1]
      temp_tomorrow <- temp_y[possible_days + 1]
      date_tomorrow <- date_y[possible_days + 1]

      # Anomalies
      cur_sim_prcp_anom <- if (!is.na(sim_prcp[prev_idx])) sim_prcp[prev_idx] - mean_mon_PRCP[cur_month] else 0
      cur_sim_temp_anom <- if (!is.na(sim_temp[prev_idx])) sim_temp[prev_idx] - mean_mon_TEMP[cur_month] else 0
      prcp_today_anom <- prcp_today - mean_mon_PRCP[cur_month]
      temp_today_anom <- temp_today - mean_mon_TEMP[cur_month]
      weights <- c(100 / sd_mon_PRCP[cur_month], 10 / sd_mon_TEMP[cur_month])

      result <- knn_sample(
        candidates = cbind(prcp_today_anom, temp_today_anom),
        target = c(cur_sim_prcp_anom, cur_sim_temp_anom),
        k = k_day,
        n = 1,
        prob = TRUE,
        weights = weights,
        seed = seed + k1 * sim_idx
      )
      result_idx <- get_result_index(result, prcp_tomorrow)

      # Save simulation values
      sim_prcp[sim_idx] <- if (!is.na(result_idx)) prcp_tomorrow[result_idx] else sim_prcp[prev_idx]
      sim_temp[sim_idx] <- if (!is.na(result_idx)) temp_tomorrow[result_idx] else sim_temp[prev_idx]
      sim_date[sim_idx] <- if (!is.na(result_idx)) date_tomorrow[result_idx] else sim_date[prev_idx]

      } # End daily loop for this year
  } # End year loop

  class(sim_date) <- "Date"
  sim_date
}
