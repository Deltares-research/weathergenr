#' K-Nearest Neighbor (KNN) Sampling from Candidates
#'
#' @description
#' Given a set of candidate vectors and a target vector, selects indices of the
#' \code{k} nearest candidates (using weighted Euclidean distance) to the target,
#' then samples \code{n} indices from these neighbors, either uniformly or with
#' rank-based or distance-based probabilities. Optionally, a random seed can be
#' set for reproducibility.
#'
#' @param candidates Numeric matrix or data frame. Each row is a candidate vector.
#' @param target Numeric vector. The target vector for comparison.
#' @param k Integer. Number of nearest neighbors to consider (k <= number of candidates).
#' @param n Integer. Number of samples to draw (default = 1).
#' @param prob Logical. If TRUE, sampling probabilities favor closer neighbors;
#'   if FALSE (default), sampling is uniform among k nearest neighbors.
#' @param weights Optional numeric vector of length equal to ncol(candidates).
#'   Feature weights for distance calculation. Default is equal weights.
#' @param seed Optional integer. If provided, sets random seed for reproducible
#'   sampling. The seed state is restored on exit.
#' @param sampling Character. Sampling method when prob = TRUE. Either "rank"
#'   (default, probability decreases with neighbor rank: 1, 1/2, 1/3, ...) or
#'   "distance" (probability based on Gaussian kernel of distances).
#'   Ignored when prob = FALSE.
#' @param bandwidth Numeric. Bandwidth parameter for distance-based sampling
#'   kernel. If NULL (default), uses median nearest neighbor distance.
#'   Only used when sampling = "distance" and prob = TRUE.
#' @param epsilon Numeric. Small constant added to distance-based probabilities
#'   to prevent zero probabilities (default = 1e-8). Only used when
#'   sampling = "distance" and prob = TRUE.
#'
#' @return
#' Integer vector of length \code{n}, giving indices of sampled candidates
#' (rows of \code{candidates}).
#'
#' @details
#' The function computes the weighted Euclidean distance between each candidate
#' and the target. The \code{k} nearest neighbors (smallest distances) are
#' identified.
#'
#' Sampling modes:
#' \describe{
#'   \item{prob = FALSE}{All k neighbors sampled with equal probability (uniform)}
#'   \item{prob = TRUE, sampling = "rank"}{Probability = 1/rank, normalized.
#'         Closest neighbor has highest probability.}
#'   \item{prob = TRUE, sampling = "distance"}{Probability based on Gaussian
#'         kernel: exp(-distance^2 / (2 * bandwidth^2))}
#' }
#'
#' Sampling is with replacement. If \code{seed} is set, the random seed will be
#' temporarily changed and restored on exit.
#'
#' @examples
#' set.seed(42)
#' candidates <- matrix(rnorm(50), ncol = 2)
#' target <- c(0, 0)
#'
#' # Sample 1 index from 5 nearest neighbors, uniform probability
#' knn_sample(candidates, target, k = 5, n = 1)
#'
#' # Sample 3 indices from 5 nearest neighbors, rank-weighted probability
#' knn_sample(candidates, target, k = 5, n = 3, prob = TRUE, seed = 123)
#'
#' # Using feature weights (weight first dimension more heavily)
#' knn_sample(candidates, target, k = 5, n = 2, weights = c(2, 1), seed = 10)
#'
#' # Distance-based sampling with custom bandwidth
#' knn_sample(candidates, target, k = 10, n = 5, prob = TRUE,
#'            sampling = "distance", bandwidth = 1.5)
#' @export
knn_sample <- function(
    candidates,
    target,
    k,
    n = 1,
    prob = FALSE,
    weights = NULL,
    seed = NULL,
    sampling = c("rank", "distance"),
    bandwidth = NULL,
    epsilon = 1e-8
) {

  sampling <- match.arg(sampling)

  # -------------------------------------------------
  # RNG handling
  # -------------------------------------------------
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }
    on.exit({ if (has_seed) .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed)
  }

  candidates <- as.matrix(candidates)
  nc <- nrow(candidates)
  p <- ncol(candidates)
  if (nc == 0) {
    stop("No candidates provided to knn_sample")
  }

  if (is.null(weights)) {
    weights <- rep(1, p)
  } else {
    if (length(weights) != p) {
      stop("Length of weights must equal number of columns in candidates.")
    }
  }

  # -------------------------------------------------
  # Weighted squared Euclidean distances
  # -------------------------------------------------
  if (p == 2) {
    d1 <- candidates[, 1] - target[1]
    d2 <- candidates[, 2] - target[2]
    d2_sq <- weights[1] * d1 * d1 + weights[2] * d2 * d2

  } else if (p <= 5) {
    diffs <- sweep(candidates, 2, target, "-")
    weighted_sq <- sweep(diffs^2, 2, weights, "*")
    d2_sq <- rowSums(weighted_sq)

  } else {
    d2_sq <- numeric(nc)
    for (j in seq_len(p)) {
      dj <- candidates[, j] - target[j]
      d2_sq <- d2_sq + weights[j] * dj * dj
    }
  }

  # -------------------------------------------------
  # Partial sorting when k << n
  # -------------------------------------------------
  k_eff <- min(k, nc)

  if (k_eff < nc * 0.2) {
    threshold <- sort(d2_sq, partial = k_eff)[k_eff]
    candidates_idx <- which(d2_sq <= threshold)
    sorted_subset_order <- order(d2_sq[candidates_idx])
    nn_indices <- candidates_idx[sorted_subset_order[seq_len(k_eff)]]
  } else {
    nn_indices <- order(d2_sq)[seq_len(k_eff)]
  }

  nn_dists <- sqrt(d2_sq[nn_indices])

  # -------------------------------------------------
  # Sampling probabilities
  # -------------------------------------------------
  if (!prob) {
    probs <- rep(1 / k_eff, k_eff)

  } else if (sampling == "rank") {
    probs <- (1 / seq_len(k_eff))
    probs <- probs / sum(probs)

  } else if (sampling == "distance") {
    if (is.null(bandwidth)) {
      bandwidth <- median(nn_dists, na.rm = TRUE)
    }

    if (!is.finite(bandwidth) || bandwidth <= 0) {
      probs <- rep(1 / k_eff, k_eff)
    } else {
      probs <- exp(-(nn_dists^2) / (2 * bandwidth^2)) + epsilon
      probs <- probs / sum(probs)
    }
  }

  # -------------------------------------------------
  # Sample neighbors
  # -------------------------------------------------
  sampled_rel <- sample.int(k_eff, n, replace = TRUE, prob = probs)
  nn_indices[sampled_rel]
}


#' Resample Daily Weather Dates Using Annual KNN and Markov Chain Logic
#'
#' @description
#' Resamples daily precipitation and temperature sequences by combining:
#' \itemize{
#'   \item annual K-nearest-neighbor (KNN) selection of observed years,
#'   \item a three-state Markov chain for wet and dry spell persistence,
#'   \item daily KNN resampling on precipitation and temperature anomalies.
#' }
#'
#' The function supports both calendar-year and water-year simulations.
#' The regime is inferred automatically from `year_start_month`:
#' if `year_start_month == 1`, calendar-year logic is used; otherwise, water-year
#' logic is assumed.
#'
#' @param sim_annual_prcp Numeric vector of length `n_years`. Synthetic
#'   annual precipitation totals generated by the annual model (e.g., WARM),
#'   indexed by simulated year.
#' @param obs_annual_prcp Numeric vector. Observed annual precipitation totals
#'   corresponding to historical water years.
#' @param obs_daily_prcp Numeric vector. Observed daily precipitation values
#'   (no leap days), aligned with `obs_dates_df`.
#' @param obs_daily_temp Numeric vector. Observed daily temperature values
#'   (no leap days), aligned with `obs_dates_df`.
#' @param year_start Integer. First simulation year (calendar year if
#'   `year_start_month == 1`, otherwise first water year).
#' @param realization_idx Integer. Realization index, used to perturb the random
#'   seed so that multiple realizations are independent.
#' @param n_years Integer. Number of simulated years.
#' @param obs_dates_df Data frame containing observed date information. Must
#'   include columns `date`, `month`, `day`, and `wyear`.
#' @param sim_dates_df Data frame containing simulated date information
#'   (no leap days). Must include columns `month`, `day`, and `wyear`.
#' @param year_start_month Integer in 1:12. First month of the simulation year.
#'   Use 1 for calendar-year simulations, or another month for water-year
#'   simulations.
#' @param annual_knn_n Integer. Number of historical years sampled in the
#'   annual KNN step.
#' @param wet_q Numeric between 0 and 1. Quantile used to define the wet
#'   threshold for daily precipitation states.
#' @param extreme_q Numeric between 0 and 1. Quantile used to define the very
#'   wet (extreme) precipitation threshold.
#' @param dry_spell_factor Numeric vector of length 12. Monthly adjustment
#'   factors controlling dry spell persistence in the Markov chain.
#' @param wet_spell_factor Numeric vector of length 12. Monthly adjustment
#'   factors controlling wet spell persistence in the Markov chain.
#' @param seed Optional integer. Base random seed for reproducibility.
#'
#' @details
#' For each simulated year, a subset of observed water years is selected
#' using annual KNN matching against `sim_annual_prcp`. Daily weather is then
#' generated sequentially:
#'
#' \enumerate{
#'   \item The first day of each simulated year is sampled from observed
#'         days matching the simulated month and day.
#'   \item Subsequent days are generated using a Markov chain to simulate
#'         wet and dry states.
#'   \item Conditional on the simulated state, candidate observed days are
#'         selected using calendar constraints and expanded search windows.
#'   \item A daily KNN step on precipitation and temperature anomalies is
#'         used to select the next-day weather values.
#' }
#'
#' In calendar-year mode (`year_start_month == 1`), transitions that cross
#' December to January are explicitly excluded to avoid cross-year
#' contamination.
#'
#' @return
#' A `Date` vector of length equal to `nrow(sim_dates_df)`, giving the
#' resampled observed dates corresponding to each simulated day.
#'
#' @importFrom stats runif
#' @export
resample_weather_dates <- function(
    sim_annual_prcp,
    obs_annual_prcp,
    obs_daily_prcp,
    obs_daily_temp,
    year_start,
    realization_idx,
    n_years,
    obs_dates_df,
    sim_dates_df,
    year_start_month = 1,
    annual_knn_n = 100,
    wet_q = 0.2,
    extreme_q = 0.8,
    dry_spell_factor = rep(1, 12),
    wet_spell_factor = rep(1, 12),
    seed = NULL
) {

  if (!is.numeric(year_start_month) || length(year_start_month) != 1L ||
      is.na(year_start_month) || !year_start_month %in% 1:12) {
    stop("year_start_month must be a single integer between 1 and 12", call. = FALSE)
  }

  # SET RNG
  base_seed <- if (is.null(seed)) NULL else seed + realization_idx

  if (!is.null(base_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }
    on.exit({ if (has_seed) .Random.seed <<- old_seed }, add = TRUE)
    set.seed(base_seed)
  }

  # Date logic and parameters
  use_water_year <- (year_start_month != 1L)
  year_month_order <- if (year_start_month == 1L) 1:12 else c(year_start_month:12, 1:(year_start_month - 1L))

  month_to_year_index <- integer(12)
  month_to_year_index[year_month_order] <- seq_along(year_month_order)

  obs_date  <- obs_dates_df$date
  obs_month <- obs_dates_df$month
  obs_day   <- obs_dates_df$day
  obs_wyear <- obs_dates_df$wyear
  obs_wyear_levels <- sort(unique(obs_dates_df$wyear))

  n_sim_day <- nrow(sim_dates_df)
  sim_month <- sim_dates_df$month
  sim_day   <- sim_dates_df$day
  sim_wyear <- sim_dates_df$wyear

  sim_daily_prcp <- numeric(n_sim_day)
  sim_daily_temp <- numeric(n_sim_day)
  sim_prcp_state <- integer(n_sim_day)
  sim_obs_date <- as.Date(rep(NA, n_sim_day))
  sim_obs_idx <- integer(n_sim_day)

  # Weights used in daily KNN sampling
  knn_weight_prcp <- 100
  knn_weight_temp <- 10

  # Constants
  offsets7  <- -3:3
  offsets61 <- -30:30

  rn_all <- runif(n_sim_day)
  k_annual <- ceiling(sqrt(length(obs_annual_prcp)))

  for (y in seq_len(n_years)) {

    year_day0_idx <- (y - 1L) * 365L + 1L
    if (year_day0_idx > n_sim_day) break

    # Annual KNN: select observed years
    obs_year_draw_idx <- knn_sample(
      candidates = obs_annual_prcp,
      target = sim_annual_prcp[y],
      k = k_annual,
      n = annual_knn_n,
      prob = TRUE,
      seed = if (is.null(base_seed)) NULL else base_seed + y
    )

    obs_year_draw <- obs_wyear_levels[obs_year_draw_idx]

    obs_idx_by_year <- split(seq_along(obs_wyear), obs_wyear)
    obs_subset_idx <- unlist(obs_idx_by_year[as.character(obs_year_draw)], use.names = FALSE)
    obs_subset_idx <- obs_subset_idx[!is.na(obs_subset_idx)]

    obs_prcp_sub  <- obs_daily_prcp[obs_subset_idx]
    obs_temp_sub  <- obs_daily_temp[obs_subset_idx]
    obs_date_sub  <- obs_date[obs_subset_idx]
    obs_month_sub <- obs_month[obs_subset_idx]
    obs_day_sub   <- obs_day[obs_subset_idx]
    obs_wyear_sub <- obs_wyear[obs_subset_idx]

    # Month-day lookup (subset-specific)
    obs_monthday_key <- paste(obs_month_sub, obs_day_sub, sep = ".")
    obs_day_lookup <- split(seq_along(obs_monthday_key), obs_monthday_key)

    # Monthly means
    obs_month_mean_prcp <- tapply(obs_prcp_sub, obs_month_sub, mean)
    obs_month_mean_temp <- tapply(obs_temp_sub, obs_month_sub, mean)
    obs_month_mean_prcp[is.na(obs_month_mean_prcp)] <- mean(obs_prcp_sub)
    obs_month_mean_temp[is.na(obs_month_mean_temp)] <- mean(obs_temp_sub)

    # Monthly standard deviations
    obs_month_sd_prcp <- tapply(obs_prcp_sub, obs_month_sub, sd)
    obs_month_sd_temp <- tapply(obs_temp_sub, obs_month_sub, sd)

    # Force full 12-month coverage
    obs_month_sd_prcp <- obs_month_sd_prcp[as.character(1:12)]
    obs_month_sd_temp <- obs_month_sd_temp[as.character(1:12)]

    # Global fallback
    obs_month_sd_prcp[is.na(obs_month_sd_prcp)] <- sd(obs_prcp_sub)
    obs_month_sd_temp[is.na(obs_month_sd_temp)] <- sd(obs_temp_sub)

    # Apply SD floor (numerical consistency)
    sd_floor_prcp <- 0.1
    sd_floor_temp <- 0.1
    obs_month_sd_prcp <- pmax(obs_month_sd_prcp, sd_floor_prcp)
    obs_month_sd_temp <- pmax(obs_month_sd_temp, sd_floor_temp)

    knn_weights_by_month <- lapply(seq_along(year_month_order), function(i) {
      c(knn_weight_prcp / obs_month_sd_prcp[year_month_order[i]],
        knn_weight_temp / obs_month_sd_temp[year_month_order[i]])
    })

    # Thresholds aligned to year_month_order
    wet_threshold_by_month_order <- sapply(year_month_order, function(m) {
      vals <- obs_prcp_sub[obs_month_sub == m]
      if (!length(vals)) NA_real_ else quantile(vals, wet_q, names = FALSE)
    })

    extreme_threshold_by_month_order <- sapply(year_month_order, function(m) {
      vals <- obs_prcp_sub[obs_month_sub == m & obs_prcp_sub > 0]
      if (!length(vals)) NA_real_ else quantile(vals, extreme_q, names = FALSE)
    })

    wet_threshold_by_month_order[is.na(wet_threshold_by_month_order)] <- quantile(obs_prcp_sub, wet_q)
    extreme_threshold_by_month_order[is.na(extreme_threshold_by_month_order)] <-
      quantile(obs_prcp_sub[obs_prcp_sub > 0], extreme_q)

    # Markov probabilities
    markov_probs <- estimate_monthly_markov_probs(
      precip_lag0 = obs_prcp_sub[-1],
      precip_lag1 = obs_prcp_sub[-length(obs_prcp_sub)],
      month_lag0  = obs_month_sub[-1],
      month_lag1  = obs_month_sub[-length(obs_prcp_sub)],
      year_lag0   = obs_wyear_sub[-1],
      year_lag1   = obs_wyear_sub[-length(obs_wyear_sub)],
      year_idx = y,
      wet_threshold = wet_threshold_by_month_order,
      extreme_threshold = extreme_threshold_by_month_order,
      month_order = year_month_order,
      sim_month = sim_month,
      sim_wyear = sim_wyear,
      sim_start_year = year_start,
      n_days_sim = n_sim_day,
      dry_spell_factor_month = dry_spell_factor,
      wet_spell_factor_month = wet_spell_factor,
      dirichlet_alpha = 1.0
    )

    # FIRST DAY OF THIS SIM YEAR
    first_month <- sim_month[year_day0_idx]
    first_month_order_idx <- match(first_month, year_month_order)
    first_day <- sim_day[year_day0_idx]

    key0 <- paste(first_month, first_day, sep = ".")
    day0_candidates <- obs_day_lookup[[key0]]

    # Fallback options
    if (!length(day0_candidates)) day0_candidates <- which(obs_month_sub == first_month)
    if (!length(day0_candidates)) day0_candidates <- seq_along(obs_prcp_sub)

    # Calendar-year safeguard: forbid observed Dec->Jan cross-year transitions
    if (!use_water_year && y > 1L) {
      prev_obs_year <- obs_wyear_sub[sim_obs_idx[year_day0_idx - 1L]]
      day0_candidates <- day0_candidates[obs_wyear_sub[day0_candidates] == prev_obs_year]
      if (!length(day0_candidates)) {
        day0_candidates <- which(obs_wyear_sub == prev_obs_year & obs_month_sub == first_month)
      }
    }

    obs_i0 <- sample(day0_candidates, 1L)
    sim_daily_prcp[year_day0_idx] <- obs_prcp_sub[obs_i0]
    sim_daily_temp[year_day0_idx] <- obs_temp_sub[obs_i0]
    sim_obs_date[year_day0_idx] <- obs_date_sub[obs_i0]
    sim_obs_idx[year_day0_idx] <- obs_i0

    # Set first day's markov state
    if (sim_daily_prcp[year_day0_idx] <= wet_threshold_by_month_order[first_month_order_idx]) {
      sim_prcp_state[year_day0_idx] <- 0L
    } else if (sim_daily_prcp[year_day0_idx] <= extreme_threshold_by_month_order[first_month_order_idx]) {
      sim_prcp_state[year_day0_idx] <- 1L
    } else {
      sim_prcp_state[year_day0_idx] <- 2L
    }

    # DAILY LOOP
    t_sim <- year_day0_idx

    for (j in 2:365) {

      t_sim <- t_sim + 1L
      if (t_sim > n_sim_day) break

      t_prev <- t_sim - 1L

      sim_prcp_state[t_sim] <- markov_next_state(
        state_prev = sim_prcp_state[t_prev],
        u_rand = rn_all[t_prev],
        idx = t_prev,
        p00 = markov_probs$p00_final,
        p01 = markov_probs$p01_final,
        p10 = markov_probs$p10_final,
        p11 = markov_probs$p11_final,
        p20 = markov_probs$p20_final,
        p21 = markov_probs$p21_final
      )

      cur_month <- sim_month[t_sim]
      cur_day   <- sim_day[t_sim]
      m_idx <- month_to_year_index[cur_month]
      if (is.na(m_idx)) m_idx <- 1L

      key <- paste(cur_month, cur_day, sep = ".")
      obs_day_candidates <- obs_day_lookup[[key]]
      if (!length(obs_day_candidates)) obs_day_candidates <- which(obs_month_sub == cur_month)

      if (!length(obs_day_candidates)) {
        sim_daily_prcp[t_sim] <- sim_daily_prcp[t_prev]
        sim_daily_temp[t_sim] <- sim_daily_temp[t_prev]
        sim_obs_date[t_sim] <- sim_obs_date[t_prev]
        next
      }

      obs_window_idx <- expand_indices(obs_day_candidates, offsets7, length(obs_prcp_sub))

      if (!length(obs_window_idx)) {
        i <- sample(obs_day_candidates, 1L)
        sim_daily_prcp[t_sim] <- obs_prcp_sub[i]
        sim_daily_temp[t_sim] <- obs_temp_sub[i]
        sim_obs_date[t_sim] <- obs_date_sub[i]
        sim_obs_idx[t_sim] <- i
        next
      }

      transition_match_pos <- match_transition_positions(
        state_from = sim_prcp_state[t_prev],
        state_to = sim_prcp_state[t_sim],
        prcp_vec = obs_prcp_sub,
        day0_idx = obs_window_idx,
        wet_threshold = wet_threshold_by_month_order[m_idx],
        extreme_threshold = extreme_threshold_by_month_order[m_idx]
      )

      if (!length(transition_match_pos)) {
        obs_window_idx <- expand_indices(obs_day_candidates, offsets61, length(obs_prcp_sub))
        transition_match_pos <- match_transition_positions(
          state_from = sim_prcp_state[t_prev],
          state_to = sim_prcp_state[t_sim],
          prcp_vec = obs_prcp_sub,
          day0_idx = obs_window_idx,
          wet_threshold = wet_threshold_by_month_order[m_idx],
          extreme_threshold = extreme_threshold_by_month_order[m_idx]
        )
      }

      if (!length(transition_match_pos)) {
        fb <- which(obs_month_sub == cur_month & (seq_along(obs_prcp_sub) + 1L) <= length(obs_prcp_sub))
        if (!use_water_year) fb <- fb[obs_wyear_sub[fb] == obs_wyear_sub[fb + 1L]]
        if (!length(fb)) {
          fb <- which((seq_along(obs_prcp_sub) + 1L) <= length(obs_prcp_sub))
          if (!use_water_year) fb <- fb[obs_wyear_sub[fb] == obs_wyear_sub[fb + 1L]]
        }

        i <- sample(fb, 1L)
        sim_daily_prcp[t_sim] <- obs_prcp_sub[i]
        sim_daily_temp[t_sim] <- obs_temp_sub[i]
        sim_obs_date[t_sim] <- obs_date_sub[i]
        sim_obs_idx[t_sim] <- i
        next
      }

      obs_day0_idx <- obs_window_idx[transition_match_pos]

      if (!use_water_year) {
        obs_day0_idx <- obs_day0_idx[obs_wyear_sub[obs_day0_idx] == obs_wyear_sub[obs_day0_idx + 1L]]
      }

      if (!length(obs_day0_idx)) {
        fb <- which(obs_month_sub == cur_month & (seq_along(obs_prcp_sub) + 1L) <= length(obs_prcp_sub))
        if (!use_water_year) fb <- fb[obs_wyear_sub[fb] == obs_wyear_sub[fb + 1L]]
        if (!length(fb)) {
          fb <- which((seq_along(obs_prcp_sub) + 1L) <= length(obs_prcp_sub))
          if (!use_water_year) fb <- fb[obs_wyear_sub[fb] == obs_wyear_sub[fb + 1L]]
        }

        i <- sample(fb, 1L)
        sim_daily_prcp[t_sim] <- obs_prcp_sub[i]
        sim_daily_temp[t_sim] <- obs_temp_sub[i]
        sim_obs_date[t_sim] <- obs_date_sub[i]
        sim_obs_idx[t_sim] <- i
        next
      }

      obs_day1_idx <- obs_day0_idx + 1L
      obs_prcp_day1 <- obs_prcp_sub[obs_day1_idx]
      obs_temp_day1 <- obs_temp_sub[obs_day1_idx]
      date_day1 <- obs_date_sub[obs_day1_idx]

      cur_sim_daily_prcp_anom <- sim_daily_prcp[t_prev] - obs_month_mean_prcp[cur_month]
      cur_sim_daily_temp_anom <- sim_daily_temp[t_prev] - obs_month_mean_temp[cur_month]

      prcp_day0_anom <- obs_prcp_sub[obs_day0_idx] - obs_month_mean_prcp[cur_month]
      temp_day0_anom <- obs_temp_sub[obs_day0_idx] - obs_month_mean_temp[cur_month]

      knn_daily_k <- max(1L, round(sqrt(length(obs_day0_idx))))
      knn_daily_k <- min(knn_daily_k, length(obs_day0_idx))

      knn_draw_pos <- knn_sample(
        candidates = cbind(prcp_day0_anom, temp_day0_anom),
        target = c(cur_sim_daily_prcp_anom, cur_sim_daily_temp_anom),
        k = knn_daily_k,
        n = 1,
        prob = TRUE,
        weights = knn_weights_by_month[[m_idx]],
        seed = if (is.null(base_seed)) NULL else base_seed + t_sim
      )

      draw_idx <- get_result_index(idx = knn_draw_pos, candidate_prcp = obs_prcp_day1)

      # Enforce calendar-year constraint AFTER KNN
      if (!use_water_year) {
        prev_obs_year <- obs_wyear_sub[sim_obs_idx[t_prev]]
        valid <- which(obs_wyear_sub[obs_day1_idx] == prev_obs_year)
        if (length(valid)) draw_idx <- draw_idx[draw_idx %in% valid]
        if (!length(draw_idx)) draw_idx <- sample(valid, 1L)
      }

      if (is.na(draw_idx) || draw_idx < 1L || draw_idx > length(obs_prcp_day1)) {
        draw_idx <- sample(seq_along(obs_prcp_day1), 1L)
      }

      sim_daily_prcp[t_sim] <- obs_prcp_day1[draw_idx]
      sim_daily_temp[t_sim] <- obs_temp_day1[draw_idx]
      sim_obs_date[t_sim] <- date_day1[draw_idx]
      sim_obs_idx[t_sim] <- obs_day1_idx[draw_idx]
    }
  }

  class(sim_obs_date) <- "Date"
  sim_obs_date
}


#' Expand Index Positions by Fixed Offsets
#'
#' @description
#' Expands a vector of base index positions by adding a set of integer offsets,
#' returning all valid index positions within a specified upper bound. This is
#' typically used to generate moving windows or neighborhood indices around
#' reference positions (e.g., expanding event days to adjacent days).
#'
#' @param base_idx Integer vector of base index positions to be expanded.
#' @param offset_vec Integer vector of offsets to apply to each element of
#'   \code{base_idx}. Offsets may be positive, negative, or zero.
#' @param n_max Integer. Maximum allowable index length; expanded indices are
#'   constrained to the range \code{1:(n_max - 1)}.
#'
#' @details
#' For each element in \code{base_idx}, all values in \code{offset_vec} are added.
#' Resulting indices that are less than or equal to zero, or for which
#' \code{index + 1 > n_max}, are discarded. The function does not enforce
#' uniqueness or sorting of the output.
#'
#' @return
#' Integer vector of expanded index positions that satisfy the validity
#' constraints.
#'
#' @examples
#' base_idx <- c(5, 10)
#' offset_vec <- -1:1
#' weathergenr:::expand_indices(base_idx, offset_vec, n_max = 15)
#'
#' @keywords internal
expand_indices <- function(base_idx, offset_vec, n_max) {
  expanded_idx <- base_idx[rep(seq_along(base_idx), each = length(offset_vec))] +
    rep(offset_vec, times = length(base_idx))
  expanded_idx[expanded_idx > 0 & (expanded_idx + 1L) <= n_max]
}


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
#' The regime is inferred from month_order: if month_order[1] == 1,
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
#' @param precip_lag0 Numeric vector. Observed daily precipitation at lag 0
#'   (current day) used to estimate state transitions.
#' @param precip_lag1 Numeric vector. Observed daily precipitation at lag 1
#'   (previous day) used to estimate state transitions.
#' @param month_lag0 Integer vector (1-12). Calendar month corresponding to
#'   precip_lag0.
#' @param month_lag1 Integer vector (1-12). Calendar month corresponding to
#'   precip_lag1.
#' @param year_lag0 Optional integer vector. Calendar year or water year
#'   corresponding to precip_lag0. If provided together with year_lag1,
#'   transitions are restricted to same-year pairs.
#' @param year_lag1 Optional integer vector. Calendar year or water year
#'   corresponding to precip_lag1.
#' @param wet_threshold Numeric vector of length 12. Monthly precipitation
#'   thresholds separating dry and wet states, aligned to month_order.
#' @param extreme_threshold Numeric vector of length 12. Monthly precipitation
#'   thresholds separating wet and very wet states, aligned to month_order.
#' @param month_order Integer vector of length 12 defining the simulated ordering
#'   of months (for example 1:12 for calendar years or c(10:12, 1:9) for
#'   October-start water years).
#' @param sim_month Integer vector. Calendar month for each simulated day.
#' @param sim_wyear Integer vector. Simulation year (calendar or water year,
#'   depending on regime) for each simulated day.
#' @param year_idx Integer. Index of the current simulated year (1-based).
#' @param sim_start_year Integer. First year of the simulation, interpreted as a
#'   calendar year or water year depending on the simulation regime.
#' @param dry_spell_factor_month Numeric vector of length 12. Monthly multiplicative
#'   factors controlling dry-state persistence.
#' @param wet_spell_factor_month Numeric vector of length 12. Monthly multiplicative
#'   factors controlling wet-state persistence.
#' @param n_days_sim Integer. Total number of simulated days. Determines the
#'   length of the returned probability vectors.
#' @param dirichlet_alpha Numeric. Base smoothing strength for transition-probability
#'   estimation. The effective smoothing applied in each month is
#'   dirichlet_alpha / sqrt(n_transitions_m), where n_transitions_m is the number of observed transitions
#'   in that month.
#'
#' @return
#' A named list of nine numeric vectors, each of length n_days_sim,
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
estimate_monthly_markov_probs <- function(
    precip_lag0,
    precip_lag1,
    month_lag0,
    month_lag1,
    year_lag0 = NULL,
    year_lag1 = NULL,
    wet_threshold,
    extreme_threshold,
    month_order,
    sim_month,
    sim_wyear,
    year_idx,
    sim_start_year,
    dry_spell_factor_month,
    wet_spell_factor_month,
    n_days_sim,
    dirichlet_alpha = 1.0
) {

  if (!is.finite(dirichlet_alpha) || dirichlet_alpha < 0) {
    stop("dirichlet_alpha must be a non-negative finite number")
  }

  use_water_year <- (month_order[1] != 1)
  sim_target_year <- if (use_water_year) (sim_start_year + year_idx) else (sim_start_year + year_idx - 1)

  use_spell_adjustment <- any(abs(dry_spell_factor_month - 1) > 1e-10) ||
    any(abs(wet_spell_factor_month - 1) > 1e-10)

  # Strong lag guard
  if (!is.null(year_lag0) && !is.null(year_lag1)) {
    keep_idx <- year_lag0 == year_lag1
    precip_lag0 <- precip_lag0[keep_idx]
    precip_lag1 <- precip_lag1[keep_idx]
    month_lag0  <- month_lag0[keep_idx]
    month_lag1  <- month_lag1[keep_idx]
  } else if (!use_water_year) {
    keep_idx <- month_lag0 >= month_lag1
    precip_lag0 <- precip_lag0[keep_idx]
    precip_lag1 <- precip_lag1[keep_idx]
    month_lag0  <- month_lag0[keep_idx]
    month_lag1  <- month_lag1[keep_idx]
  }

  # Map each observation to month index
  month_order_idx_lag1 <- match(month_lag1, month_order)
  month_order_idx_lag0 <- match(month_lag0, month_order)

  wet_threshold_lag1 <- wet_threshold[month_order_idx_lag1]
  extreme_threshold_lag1 <- extreme_threshold[month_order_idx_lag1]
  wet_threshold_lag0 <- wet_threshold[month_order_idx_lag0]
  extreme_threshold_lag0 <- extreme_threshold[month_order_idx_lag0]

  # Handle non-finite thresholds (fallback to global)
  if (any(!is.finite(wet_threshold_lag1)) || any(!is.finite(extreme_threshold_lag1))) {
    wet_threshold_global <- quantile(precip_lag1, 0.2, na.rm = TRUE)
    pos <- precip_lag1[precip_lag1 > 0]
    extreme_threshold_global <- quantile(if (length(pos)) pos else precip_lag1, 0.8, na.rm = TRUE)

    wet_threshold_lag1[!is.finite(wet_threshold_lag1)] <- wet_threshold_global
    extreme_threshold_lag1[!is.finite(extreme_threshold_lag1)] <- extreme_threshold_global
    wet_threshold_lag0[!is.finite(wet_threshold_lag0)] <- wet_threshold_global
    extreme_threshold_lag0[!is.finite(extreme_threshold_lag0)] <- extreme_threshold_global
  }

  state_lag1 <- ifelse(precip_lag1 <= wet_threshold_lag1, 0L,
                       ifelse(precip_lag1 <= extreme_threshold_lag1, 1L, 2L))

  state_lag0 <- ifelse(precip_lag0 <= wet_threshold_lag0, 0L,
                       ifelse(precip_lag0 <= extreme_threshold_lag0, 1L, 2L))

  p00_final <- p01_final <- p02_final <- rep(NA_real_, n_days_sim)
  p10_final <- p11_final <- p12_final <- rep(NA_real_, n_days_sim)
  p20_final <- p21_final <- p22_final <- rep(NA_real_, n_days_sim)

  for (m in seq_along(month_order)) {

    mm <- month_order[m]
    r <- which(sim_month == mm & sim_wyear == sim_target_year)
    if (!length(r)) next

    x <- which(month_lag1 == mm)
    n_transitions_m <- length(x)
    dirichlet_alpha_m_eff <- if (n_transitions_m > 0L) dirichlet_alpha / sqrt(n_transitions_m) else dirichlet_alpha

    if (!length(x)) {
      p00_final[r] <- 1; p01_final[r] <- 0; p02_final[r] <- 0
      p10_final[r] <- 1; p11_final[r] <- 0; p12_final[r] <- 0
      p20_final[r] <- 1; p21_final[r] <- 0; p22_final[r] <- 0
      next
    }

    s1 <- state_lag1[x]
    s0 <- state_lag0[x]

    n00 <- sum(s1 == 0 & s0 == 0)
    n01 <- sum(s1 == 0 & s0 == 1)
    n02 <- sum(s1 == 0 & s0 == 2)

    n10 <- sum(s1 == 1 & s0 == 0)
    n11 <- sum(s1 == 1 & s0 == 1)
    n12 <- sum(s1 == 1 & s0 == 2)

    n20 <- sum(s1 == 2 & s0 == 0)
    n21 <- sum(s1 == 2 & s0 == 1)
    n22 <- sum(s1 == 2 & s0 == 2)

    p_dry <- c(n00, n01, n02) + dirichlet_alpha_m_eff
    p_wet <- c(n10, n11, n12) + dirichlet_alpha_m_eff
    p_vwt <- c(n20, n21, n22) + dirichlet_alpha_m_eff

    p_dry <- p_dry / sum(p_dry)
    p_wet <- p_wet / sum(p_wet)
    p_vwt <- p_vwt / sum(p_vwt)

    if (use_spell_adjustment) {
      dry_factor_m <- dry_spell_factor_month[m]
      wet_factor_m <- wet_spell_factor_month[m]

      if (!is.finite(dry_factor_m) || dry_factor_m <= 0) dry_factor_m <- 1
      if (!is.finite(wet_factor_m) || wet_factor_m <= 0) wet_factor_m <- 1

      if (abs(dry_factor_m - 1) > 1e-10) {
        p_dry <- normalize_probs(
          prob = c(p_dry[1], p_dry[2] / dry_factor_m, p_dry[3] / dry_factor_m),
          fallback_prob = c(1, 0, 0)
        )
      }

      if (abs(wet_factor_m - 1) > 1e-10) {
        p_wet <- normalize_probs(
          prob = c(p_wet[1] / wet_factor_m, p_wet[2], p_wet[3]),
          fallback_prob = c(0, 1, 0)
        )
      }
    }

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


#' Normalize a Probability Vector with Robust Fallback Handling
#'
#' @description
#' Normalizes a numeric vector to sum to one, treating non-finite and negative
#' values as zero. If the total mass is zero after cleaning, a fallback
#' distribution is returned instead.
#'
#' @param prob Numeric vector of (unnormalized) probabilities or weights.
#'   Non-finite values (\code{NA}, \code{NaN}, \code{Inf}) and negative values
#'   are set to zero prior to normalization.
#' @param fallback_prob Optional numeric vector to return when \code{prob} has zero
#'   total mass after cleaning. If \code{NULL} (default), a uniform distribution
#'   of length \code{length(prob)} is returned.
#'
#' @details
#' This helper is intended for internal use in stochastic or probabilistic
#' workflows where degenerate or numerically unstable probability vectors may
#' arise. The function does not check that \code{fallback_prob} sums to one; this is
#' assumed to be ensured by the caller.
#'
#' @return
#' Numeric vector of probabilities summing to one. If a fallback is used, the
#' returned vector is either \code{fallback_prob} or a uniform distribution.
#'
#' @examples
#' \dontrun{
#' prob <- c(0.2, NA, -0.1, 0.5)
#' normalize_probs(prob)
#'
#' normalize_probs(c(0, 0, 0))
#'
#' normalize_probs(c(0, 0, 0), fallback_prob = c(0.7, 0.2, 0.1))
#' }
#'
#' @keywords internal
normalize_probs <- function(prob, fallback_prob = NULL) {
  prob[!is.finite(prob) | prob < 0] <- 0
  s <- sum(prob)
  if (s > 0) {
    prob / s
  } else {
    if (is.null(fallback_prob)) {
      rep(1 / length(prob), length(prob))
    } else {
      fallback_prob
    }
  }
}



#' Determine Next State in a First-Order Markov Chain
#'
#' Given the previous state and a random number `u_rand`, this function returns the next
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
#' @param state_prev Integer. The current (previous) state (0, 1, or 2).
#' @param u_rand Numeric. A random number in [0, 1] used to sample the next state.
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
#' set.seed(123)
#' u_rand <- runif(1)
#' markov_next_state(
#'   state_prev = 1,
#'   u_rand = u_rand,
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
markov_next_state <- function(state_prev, u_rand, idx, p00, p01, p10, p11, p20, p21) {

  if (!is.finite(u_rand)) u_rand <- 0
  if (u_rand < 0) u_rand <- 0
  if (u_rand > 1) u_rand <- 1

  if (!(state_prev %in% c(0L, 1L, 2L))) {
    state_prev <- 0L
  }

  n <- length(p00)
  if (idx < 1L) idx <- 1L
  if (idx > n)  idx <- n

  if (state_prev == 0L) {
    p_to0 <- p00[idx]
    p_to1 <- p01[idx]
  } else if (state_prev == 1L) {
    p_to0 <- p10[idx]
    p_to1 <- p11[idx]
  } else {
    p_to0 <- p20[idx]
    p_to1 <- p21[idx]
  }

  if (!is.finite(p_to0)) p_to0 <- 0
  if (!is.finite(p_to1)) p_to1 <- 0

  if (p_to0 < 0) p_to0 <- 0
  if (p_to1 < 0) p_to1 <- 0

  row_sum <- p_to0 + p_to1
  if (row_sum > 1) {
    if (row_sum > 1.01) {
      warning(sprintf(
        "Invalid transition probabilities: sum = %.3f at index %s.",
        row_sum,
        format(idx, big.mark = ",")
      ))
    }
    p_to0 <- p_to0 / row_sum
    p_to1 <- p_to1 / row_sum
  }

  if (u_rand < p_to0) {
    0L
  } else if (u_rand < (p_to0 + p_to1)) {
    1L
  } else {
    2L
  }
}


#' Get Positions of a Specific Occurrence State Transition Within Candidate Indices
#'
#' Returns the positions within `day0_idx` where a specified state transition occurs,
#' based on wet/dry/extreme day thresholds.
#'
#' Occurrence states are defined as:
#' - 0: Dry day (precipitation <= wet_threshold)
#' - 1: Wet day (precipitation > wet_threshold and <= extreme_threshold)
#' - 2: Extreme wet day (precipitation > extreme_threshold)
#'
#' @param state_from Integer. The current occurrence state (0 = dry, 1 = wet, 2 = extreme).
#' @param state_to Integer. The next occurrence state (0 = dry, 1 = wet, 2 = extreme).
#' @param prcp_vec Numeric vector of precipitation values.
#' @param day0_idx Integer vector of indices representing candidate "day 0" positions
#'   in the time series.
#' @param wet_threshold Numeric. Threshold separating dry and wet days.
#' @param extreme_threshold Numeric. Threshold above which days are considered extreme.
#'
#' @return
#' Integer vector of positions within `day0_idx` where the transition from
#' `state_from` to `state_to` occurs.
#'
#' @examples
#' prcp_vec <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
#' day0_idx <- 1:(length(prcp_vec) - 1)
#' wet_threshold <- 1
#' extreme_threshold <- 20
#'
#' match_transition_positions(0, 1, prcp_vec, day0_idx, wet_threshold, extreme_threshold)
#' match_transition_positions(1, 2, prcp_vec, day0_idx, wet_threshold, extreme_threshold)
#' match_transition_positions(2, 0, prcp_vec, day0_idx, wet_threshold, extreme_threshold)
#'
#' @export
match_transition_positions <- function(
    state_from,
    state_to,
    prcp_vec,
    day0_idx,
    wet_threshold,
    extreme_threshold
) {

  prcp_day0 <- prcp_vec[day0_idx]
  prcp_day1 <- prcp_vec[day0_idx + 1L]

  state_day0 <- ifelse(prcp_day0 <= wet_threshold, 0L,
                       ifelse(prcp_day0 <= extreme_threshold, 1L, 2L))

  state_day1 <- ifelse(prcp_day1 <= wet_threshold, 0L,
                       ifelse(prcp_day1 <= extreme_threshold, 1L, 2L))

  which(state_day0 == state_from & state_day1 == state_to)
}


#' Safely Get or Sample an Index from a Vector
#'
#' Returns a valid index from `candidate_prcp` based on the provided `idx`.
#' If `idx` is missing, out of bounds, or invalid, a random index is sampled instead.
#' If `candidate_prcp` is empty, `NA_integer_` is returned.
#'
#' @param idx Integer. Proposed index for retrieving an element from `candidate_prcp`.
#' @param candidate_prcp Numeric vector. A set of candidate precipitation values (e.g., for the next day).
#'
#' @return Integer. A valid index from `candidate_prcp`, or `NA_integer_` if `candidate_prcp` is empty.
#'
#' @examples
#' candidate_prcp <- c(0, 5, 10, 20)
#' get_result_index(2, candidate_prcp)
#'
#' set.seed(1)
#' get_result_index(10, candidate_prcp)
#'
#' get_result_index(NA, candidate_prcp)
#' get_result_index(1, numeric(0))
#'
#' @export
get_result_index <- function(idx, candidate_prcp) {
  if (is.na(idx) || length(idx) == 0 || idx < 1 || idx > length(candidate_prcp)) {
    if (length(candidate_prcp) > 0) {
      sample(seq_along(candidate_prcp), 1)
    } else {
      NA_integer_
    }
  } else {
    idx
  }
}
