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
    START_YEAR_SIM,
    k1,
    ymax,
    dates.d,
    sim.dates.d,
    YEAR_D,
    month.start = 1,
    knn.annual.sample.num = 50,
    wet.quantile = 0.2,
    extreme.quantile = 0.8,
    dry.spell.change = rep(1, 12),
    wet.spell.change = rep(1, 12),
    seed = NULL) {
  ######################## --- Input Validation --- ############################

  # stopifnot(
  #   !is.null(PRCP_FINAL_ANNUAL_SIM),
  #   !is.null(ANNUAL_PRCP),
  #   !is.null(PRCP),
  #   !is.null(TEMP),
  #   !is.null(START_YEAR_SIM),
  #   !is.null(dates.d),
  #   !is.null(sim.dates.d),
  #   !is.null(YEAR_D),
  #   !is.null(month.start)
  # )
  #
  # # Numeric checks
  # if (!is.numeric(PRCP_FINAL_ANNUAL_SIM) || !is.numeric(ANNUAL_PRCP))
  #   stop("PRCP_FINAL_ANNUAL_SIM and ANNUAL_PRCP must be numeric vectors.")
  # if (!is.numeric(PRCP) || !is.numeric(TEMP))
  #   stop("PRCP and TEMP must be numeric vectors.")
  # if (!is.numeric(wet.quantile) || wet.quantile < 0 || wet.quantile > 1)
  #   stop("wet.quantile must be between 0 and 1.")
  # if (!is.numeric(extreme.quantile) || extreme.quantile < 0 || extreme.quantile > 1)
  #   stop("extreme.quantile must be between 0 and 1.")
  # if (!is.numeric(knn.annual.sample.num) || knn.annual.sample.num < 1)
  #   stop("knn.annual.sample.num must be a positive integer.")
  # if (!is.numeric(seed) && !is.null(seed))
  #   stop("seed must be numeric or NULL.")
  #
  # # Data frame checks
  # if (!is.data.frame(dates.d) || !is.data.frame(sim.dates.d))
  #   stop("dates.d and sim.dates.d must be data.frames.")
  #
  # # Check required columns
  # required_cols <- c("date", "month", "day", "wyear")
  # missing_cols_dates <- setdiff(required_cols, names(dates.d))
  # missing_cols_sim   <- setdiff(c("month", "day", "wyear"), names(sim.dates.d))
  # if (length(missing_cols_dates) > 0)
  #   stop(paste("dates.d is missing columns:", paste(missing_cols_dates, collapse = ", ")))
  # if (length(missing_cols_sim) > 0)
  #   stop(paste("sim.dates.d is missing columns:", paste(missing_cols_sim, collapse = ", ")))
  #
  # # Check lengths (for main time series vectors)
  # expected_length <- nrow(dates.d)
  # if (length(PRCP) != expected_length)
  #   stop("Length of PRCP does not match dates.d.")
  # if (length(TEMP) != expected_length)
  #   stop("Length of TEMP does not match dates.d.")
  # if (length(YEAR_D) != expected_length)
  #   stop("Length of YEAR_D does not match dates.d.")

  ######################## --- LOCAL RANDOM SEED --- ##########################

  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit(
      {
        .Random.seed <<- old_seed
      },
      add = TRUE
    )
    set.seed(seed + k1)
  }

  # Workaround for rlang warning
  month <- day <- wyear <- 0

  ####################################################################################################

  # Prepare month list and reference vectors
  month_list <- if (month.start == 1) 1:12 else c(month.start:12, 1:(month.start - 1))
  WATER_YEAR_A <- dplyr::filter(dates.d, month == month.start & day == 1) %>% dplyr::pull(wyear)
  DATE_D <- dates.d$date
  MONTH_D <- dates.d$month
  WATER_YEAR_D <- dates.d$wyear
  MONTH_DAY_D <- dates.d[, c("month", "day")]
  MONTH_SIM <- sim.dates.d$month
  DAY_SIM <- sim.dates.d$day
  WATER_YEAR_SIM <- sim.dates.d$wyear
  SIM_LENGTH <- length(MONTH_SIM)
  water_year_start <- dates.d$wyear[1]

  # Initialize Markov transition arrays and simulated states
  p00_final <- array(NA, SIM_LENGTH)
  p01_final <- array(NA, SIM_LENGTH)
  p02_final <- array(NA, SIM_LENGTH)
  p10_final <- array(NA, SIM_LENGTH)
  p11_final <- array(NA, SIM_LENGTH)
  p12_final <- array(NA, SIM_LENGTH)
  p20_final <- array(NA, SIM_LENGTH)
  p21_final <- array(NA, SIM_LENGTH)
  p22_final <- array(NA, SIM_LENGTH)
  OCCURENCES <- array(0, c(SIM_LENGTH))
  SIM_PRCP <- array(0, c(SIM_LENGTH))
  SIM_TEMP <- array(25, c(SIM_LENGTH))

  SIM_DATE <- array(as.Date(paste(water_year_start + 1, month.start, "01", sep = "-")), SIM_LENGTH)
  rn_all <- stats::runif(SIM_LENGTH, 0, 1)
  kk <- max(round(sqrt(length(ANNUAL_PRCP)), 0), round(length(ANNUAL_PRCP), 0) * .5)
  count <- 1

  # ----- 3. Simulate years ----------------------------------------------------

  # Global mean and standard deviaton
  mean_monthly_TEMP <- sapply(1:12, function(x) mean(TEMP[which(MONTH_D == x)]))
  mean_monthly_PRCP <- sapply(1:12, function(x) mean(PRCP[which(MONTH_D == x)]))

  sd_monthly_TEMP <- sapply(1:12, function(x) stats::sd(TEMP[which(MONTH_D == x)]))
  sd_monthly_PRCP <- sapply(1:12, function(x) stats::sd(PRCP[which(MONTH_D == x)]))


  for (y in seq_len(ymax)) {
    # - For each simulated year, select similar years from the observed record via KNN
    # - Calculate monthly precipitation thresholds (wet, extreme)
    # - Calculate Markov transition probabilities
    # - For each simulated day, use Markov Chain to determine state, then sample a day from historical record

    # --- 3a. Sample annual years ---
    sim_annual_prcp <- PRCP_FINAL_ANNUAL_SIM[y]

    cur_year_index <- knn_sample(
      candidates = ANNUAL_PRCP,
      target = sim_annual_prcp,
      k = kk,
      n = knn.annual.sample.num,
      prob = TRUE,
      weights = NULL,
      seed = seed + k1 * y
    )

    CUR_YEARS <- WATER_YEAR_A[cur_year_index]

    # # Find indices of days in all sampled years in CUR_YEARS
    conditional_selection <- unlist(lapply(CUR_YEARS, function(x) which(WATER_YEAR_D == x)))

    # Find all variables and date indices in the conditional selection
    PRCP_CURRENT <- PRCP[conditional_selection]
    TEMP_CURRENT <- TEMP[conditional_selection]
    DATE_D_CURRENT <- DATE_D[conditional_selection]
    MONTH_D_CURRENT <- MONTH_D[conditional_selection]
    YEAR_D_CURRENT <- YEAR_D[conditional_selection]
    MONTH_DAY_D_CURRENT <- MONTH_DAY_D[conditional_selection, ]

    # --- 3b. Thresholds ---
    wet_threshold <- rep(stats::quantile(PRCP_CURRENT, wet.quantile, names = FALSE), 12)
    extreme_threshold <- sapply(1:12, function(m) {
      stats::quantile(PRCP_CURRENT[which(MONTH_D_CURRENT == month_list[m] & PRCP_CURRENT > wet_threshold[m])],
        extreme.quantile,
        names = F
      )
    })

    ############################################################################

    # Define lagged variables on daily time-series (for current year set)
    PRCP_LAG0 <- PRCP_CURRENT[2:length(PRCP_CURRENT)]
    PRCP_LAG1 <- PRCP_CURRENT[1:(length(PRCP_CURRENT) - 1)]
    MONTH_LAG0 <- MONTH_D_CURRENT[2:length(PRCP_CURRENT)]
    MONTH_LAG1 <- MONTH_D_CURRENT[1:(length(PRCP_CURRENT) - 1)]

    # --- 3d. Vectorized Markov transition probability calculation ---
    probs <- calculateMarkovProbs(
      PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
      wet_threshold, extreme_threshold, month_list,
      MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
      dry.spell.change, wet.spell.change, SIM_LENGTH
    )

    p00_final <- probs$p00_final
    p01_final <- probs$p01_final
    p02_final <- probs$p02_final
    p10_final <- probs$p10_final
    p11_final <- probs$p11_final
    p12_final <- probs$p12_final
    p20_final <- probs$p20_final
    p21_final <- probs$p21_final
    p22_final <- probs$p22_final

    ############################################################################
    ############################################################################

    # At the top of your function/script (once)
    lookup_cur_day <- split(
      seq_len(nrow(MONTH_DAY_D_CURRENT)),
      paste(MONTH_DAY_D_CURRENT$month, MONTH_DAY_D_CURRENT$day, sep = ".")
    )


    # MARKOV-CHAIN AND DAILY KNN SAMPLING.......................................
    for (j in 1:365) {
      count <- count + 1
      idx <- count - 1 # previous day index
      if (count > SIM_LENGTH) break

      # --- Markov Chain Transition ---
      OCCURENCES[count] <- markov_next_state(
        OCCURENCES[idx], rn_all[idx], idx, p00_final, p01_final, p10_final,
        p11_final, p20_final, p21_final
      )

      # Current day's month of the year and day of the year indices
      m <- MONTH_SIM[idx]
      mmm <- which(month_list == m)
      d <- DAY_SIM[idx]

      # Current and next day occurrences
      cur_OCC <- OCCURENCES[(idx)]
      next_OCC <- OCCURENCES[(count)]

      # Subset non-zero days?
      cur_day <- lookup_cur_day[[paste(m, d, sep = ".")]]
      cur_day <- c((cur_day - 3), (cur_day - 2), (cur_day - 1), cur_day, (cur_day + 1), (cur_day + 2), (cur_day + 3))
      cur_day <- subset(cur_day, cur_day > 0)

      state_idx <- get_state_indices(
        cur_OCC, next_OCC, PRCP_CURRENT, cur_day,
        wet_threshold[m], extreme_threshold[mmm]
      )

      # --- Expand window if no suitable days found ---
      if (length(state_idx) == 0) {
        cur_day <- lookup_cur_day[[paste(m, d, sep = ".")]]
        cur_day_final <- array(NA, length(cur_day) * 61)
        cur_day_window <- seq(-30, 30)
        ncur_day <- length(cur_day)

        for (cc in 1:61) {
          cur_day_final[(1 + ncur_day * (cc - 1)):(ncur_day + ncur_day * (cc - 1))] <- (cur_day + cur_day_window[cc])
        }
        cur_day <- subset(cur_day_final, cur_day_final > 0)
        state_idx <- get_state_indices(
          cur_OCC, next_OCC, PRCP_CURRENT, cur_day,
          wet_threshold[m], extreme_threshold[mmm]
        )
      }

      possible_days <- cur_day[state_idx]

      # --- KNN Sampling Preparation ---
      PRCP_TODAY <- PRCP_CURRENT[possible_days]
      TEMP_TODAY <- TEMP_CURRENT[possible_days]
      PRCP_TOMORROW <- PRCP_CURRENT[possible_days + 1]
      TEMP_TOMORROW <- TEMP_CURRENT[possible_days + 1]
      DATE_TOMORROW <- DATE_D_CURRENT[possible_days + 1]

      k <- round(sqrt(length(possible_days)))

      cur_sim_PRCP_anom <- SIM_PRCP[(idx)] - mean_monthly_PRCP[m]
      cur_sim_TEMP_anom <- SIM_TEMP[(idx)] - mean_monthly_TEMP[m]
      PRCP_TODAY_anom <- PRCP_TODAY - mean_monthly_PRCP[m]
      TEMP_TODAY_anom <- TEMP_TODAY - mean_monthly_TEMP[m]

      # Weights
      RESULT <- knn_sample(
        candidates = cbind(PRCP_TODAY_anom, TEMP_TODAY_anom),
        target = c(cur_sim_PRCP_anom, cur_sim_TEMP_anom),
        k = k,
        n = 1,
        prob = TRUE,
        weights = c(100 / sd_monthly_PRCP[m], 10 / sd_monthly_TEMP[m]),
        seed = seed + k1 * count
      )

      result_idx <- get_result_index(RESULT, PRCP_TOMORROW)

      # --- Assign simulated values ---
      SIM_PRCP[count] <- if (!is.na(result_idx)) PRCP_TOMORROW[result_idx] else NA_real_
      SIM_TEMP[count] <- if (!is.na(result_idx)) TEMP_TOMORROW[result_idx] else NA_real_
      SIM_DATE[count] <- if (!is.na(result_idx)) DATE_TOMORROW[result_idx] else as.Date(NA)
    } # daily-counter close
  } # year-counter close

  class(SIM_DATE) <- "Date"

  return(SIM_DATE)
}
