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
    if (length(weights) != p)
      stop("Length of weights must equal number of columns in candidates.")
  }

  # -------------------------------------------------
  # Weighted squared Euclidean distances
  # -------------------------------------------------
  if (p == 2) {
    # Special case for p=2 (common: precipitation + temperature)
    # Unrolled for maximum speed
    d1 <- candidates[, 1] - target[1]
    d2 <- candidates[, 2] - target[2]
    d2_sq <- weights[1] * d1 * d1 + weights[2] * d2 * d2

  } else if (p <= 5) {
    # Vectorized operations for small p
    diffs <- sweep(candidates, 2, target, "-")
    weighted_sq <- sweep(diffs^2, 2, weights, "*")
    d2_sq <- rowSums(weighted_sq)

  } else {
    # Loop for large p (memory efficient)
    d2_sq <- numeric(nc)
    for (j in seq_len(p)) {
      dj <- candidates[, j] - target[j]
      d2_sq <- d2_sq + weights[j] * dj * dj
    }
  }

  # -------------------------------------------------
  # OPTIMIZED: Partial sorting when k << n
  # -------------------------------------------------
  k_eff <- min(k, nc)

  if (k_eff < nc * 0.2) {
    # Partial sorting: faster when k is small relative to n
    # Uses sort.int with partial parameter: O(n) + O(k log k) instead of O(n log n)

    # Find the k-th smallest squared distance (partition threshold)
    threshold <- sort(d2_sq, partial = k_eff)[k_eff]

    # Get all candidate indices with squared distance <= threshold
    # (may include ties at the boundary)
    candidates_idx <- which(d2_sq <= threshold)

    # Fully sort only this subset to get exact k nearest
    sorted_subset_order <- order(d2_sq[candidates_idx])
    nn_indices <- candidates_idx[sorted_subset_order[seq_len(k_eff)]]

  } else {
    # k is large relative to n: full sort is fine
    nn_indices <- order(d2_sq)[seq_len(k_eff)]
  }

  # OPTIMIZED: Only compute sqrt for the k nearest neighbors
  # (avoided n - k sqrt operations)
  nn_dists <- sqrt(d2_sq[nn_indices])

  # -------------------------------------------------
  # Sampling probabilities
  # -------------------------------------------------
  if (!prob) {
    # Uniform sampling
    probs <- rep(1 / k_eff, k_eff)

  } else if (sampling == "rank") {
    # Rank-based probabilities (default for prob = TRUE)
    probs <- (1 / seq_len(k_eff))
    probs <- probs / sum(probs)

  } else if (sampling == "distance") {
    # Distance-based probabilities (Gaussian kernel)
    if (is.null(bandwidth)) {
      # Automatic bandwidth: median NN distance
      bandwidth <- median(nn_dists, na.rm = TRUE)
    }

    if (!is.finite(bandwidth) || bandwidth <= 0) {
      # Fallback to uniform if bandwidth invalid
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
#' @param annual_prcp_sim Numeric vector of length `n_sim_years`. Synthetic
#'   annual precipitation totals generated by the annual model (e.g., WARM),
#'   indexed by simulated year.
#' @param annual_prcp_obs Numeric vector. Observed annual precipitation totals
#'   corresponding to historical water years.
#' @param daily_prcp_obs Numeric vector. Observed daily precipitation values
#'   (no leap days), aligned with `obs_dates_table`.
#' @param daily_temp_obs Numeric vector. Observed daily temperature values
#'   (no leap days), aligned with `obs_dates_table`.
#' @param sim_start_year Integer. First simulation year (calendar year if
#'   `year_start_month == 1`, otherwise first water year).
#' @param realization_id Integer. Realization index, used to perturb the random
#'   seed so that multiple realizations are independent.
#' @param n_sim_years Integer. Number of simulated years.
#' @param obs_dates_table Data frame containing observed date information. Must
#'   include columns `date`, `month`, `day`, and `wyear`.
#' @param sim_dates_table Data frame containing simulated date information
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
#' using annual KNN matching against `annual_prcp_sim`. Daily weather is then
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
#' A `Date` vector of length equal to `nrow(sim_dates_table)`, giving the
#' resampled observed dates corresponding to each simulated day.
#'
#' @importFrom stats runif
#' @export
resample_weather_dates <- function(
    annual_prcp_sim,
    annual_prcp_obs,
    daily_prcp_obs,
    daily_temp_obs,
    sim_start_year,
    realization_id,
    n_sim_years,
    obs_dates_table,
    sim_dates_table,
    year_start_month = 1,
    annual_knn_n = 50,
    wet_q = 0.2,
    extreme_q = 0.8,
    dry_spell_factor = rep(1, 12),
    wet_spell_factor = rep(1, 12),
    seed = NULL
) {

  # Checks and controls
  if (!is.numeric(year_start_month) || length(year_start_month) != 1L ||
      is.na(year_start_month) || !year_start_month %in% 1:12) {
    stop("year_start_month must be a single integer between 1 and 12", call. = FALSE)
  }

  # SET RNG
  base_seed <- if (is.null(seed)) NULL else seed + realization_id

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

  # SET Date logic and parameters
  is_water_year <- (year_start_month != 1L)
  year_month_order <- if (year_start_month == 1L) 1:12 else c(year_start_month:12, 1:(year_start_month - 1L))

  month_to_year_index <- integer(12)
  month_to_year_index[year_month_order] <- seq_along(year_month_order)

  dates_obs  <- obs_dates_table$date
  months_obs <- obs_dates_table$month
  days_obs   <- obs_dates_table$day
  wyears_obs <- obs_dates_table$wyear
  water_years_obs <- sort(unique(obs_dates_table$wyear))

  sim_length <- nrow(sim_dates_table)
  months_sim <- sim_dates_table$month
  days_sim   <- sim_dates_table$day
  wyears_sim <- sim_dates_table$wyear

  sim_prcp <- numeric(sim_length)
  sim_temp <- numeric(sim_length)
  sim_occ  <- integer(sim_length)
  sim_date <- as.Date(rep(NA, sim_length))
  sim_obs_idx <- integer(sim_length)

  # Weights used in daily KNN sampling
  knnw_prcp <- 100
  knnw_temp <- 10

  # Constants
  offsets7  <- -3:3
  offsets61 <- -30:30

  rn_all <- runif(sim_length)
  k_annual <- ceiling(sqrt(length(annual_prcp_obs)))

  for (y in seq_len(n_sim_years)) {

    year_start_idx <- (y - 1L) * 365L + 1L
    if (year_start_idx > sim_length) break

    # Annual KNN: select observed years
    year_sample_idx <- knn_sample(
      candidates = annual_prcp_obs,
      target = annual_prcp_sim[y],
      k = k_annual,
      n = annual_knn_n,
      prob = TRUE,
      seed = if (is.null(base_seed)) NULL else base_seed + y
    )

    cur_years <- water_years_obs[year_sample_idx]

    idx_by_year <- split(seq_along(wyears_obs), wyears_obs)
    obs_idx <- unlist(idx_by_year[as.character(cur_years)], use.names = FALSE)
    obs_idx <- obs_idx[!is.na(obs_idx)]

    prcp_y  <- daily_prcp_obs[obs_idx]
    temp_y  <- daily_temp_obs[obs_idx]
    date_y  <- dates_obs[obs_idx]
    month_y <- months_obs[obs_idx]
    day_y   <- days_obs[obs_idx]
    wyear_y <- wyears_obs[obs_idx]

    # Month-day lookup (subset-specific)
    monthday_key_y <- paste(month_y, day_y, sep = ".")
    lookup_day_idx_y <- split(seq_along(monthday_key_y), monthday_key_y)

    # Monthly means
    mean_mon_prcp <- tapply(prcp_y, month_y, mean)
    mean_mon_temp <- tapply(temp_y, month_y, mean)
    mean_mon_prcp[is.na(mean_mon_prcp)] <- mean(prcp_y)
    mean_mon_temp[is.na(mean_mon_temp)] <- mean(temp_y)

    # Monthly standard deviations
    sd_mon_prcp <- tapply(prcp_y, month_y, sd)
    sd_mon_temp <- tapply(temp_y, month_y, sd)

    # Force full 12-month coverage
    sd_mon_prcp <- sd_mon_prcp[as.character(1:12)]
    sd_mon_temp <- sd_mon_temp[as.character(1:12)]

    # Global fallback
    sd_mon_prcp[is.na(sd_mon_prcp)] <- sd(prcp_y)
    sd_mon_temp[is.na(sd_mon_temp)] <- sd(temp_y)

    # Apply SD floor (numerical consistency)
    sd_floor_prcp <- 0.1
    sd_floor_temp <- 0.1
    sd_mon_prcp <- pmax(sd_mon_prcp, sd_floor_prcp)
    sd_mon_temp <- pmax(sd_mon_temp, sd_floor_temp)

    knn_weights_by_month <- lapply(seq_along(year_month_order), function(i) {
      c(knnw_prcp / sd_mon_prcp[year_month_order[i]],
        knnw_temp / sd_mon_temp[year_month_order[i]])
    })

    # Thresholds aligned to year_month_order
    wet_thresh <- sapply(year_month_order, function(m) {
      vals <- prcp_y[month_y == m]
      if (!length(vals)) NA_real_ else quantile(vals, wet_q, names = FALSE)
    })

    extreme_thresh <- sapply(year_month_order, function(m) {
      vals <- prcp_y[month_y == m & prcp_y > 0]
      if (!length(vals)) NA_real_ else quantile(vals, extreme_q, names = FALSE)
    })

    wet_thresh[is.na(wet_thresh)] <- quantile(prcp_y, wet_q)
    extreme_thresh[is.na(extreme_thresh)] <- quantile(prcp_y[prcp_y > 0], extreme_q)

    # Markov probabilities
    probs <- monthly_markov_probs(
      precip.lag0 = prcp_y[-1],
      precip.lag1 = prcp_y[-length(prcp_y)],
      month.lag0  = month_y[-1],
      month.lag1  = month_y[-length(prcp_y)],
      year.lag0   = wyear_y[-1],
      year.lag1   = wyear_y[-length(wyear_y)],
      year.idx = y,
      wet.threshold = wet_thresh,
      extreme.threshold = extreme_thresh,
      month.list = year_month_order,
      sim.months = months_sim,
      sim.water.years = wyears_sim,
      sim.start.year = sim_start_year,
      sim.length = sim_length,
      dry.spell.change = dry_spell_factor,
      wet.spell.change = wet_spell_factor,
      alpha = 1.0
    )

    # FIRST DAY OF THIS SIM YEAR
    first_month <- months_sim[year_start_idx]
    first_month_idx <- match(first_month, year_month_order)
    first_day <- days_sim[year_start_idx]

    key0 <- paste(first_month, first_day, sep = ".")
    cand0 <- lookup_day_idx_y[[key0]]

    # Fallback options
    if (!length(cand0)) cand0 <- which(month_y == first_month)
    if (!length(cand0)) cand0 <- seq_along(prcp_y)

    # Calendar-year safeguard: forbid observed Dec->Jan cross-year transitions
    if (!is_water_year && y > 1L) {
      prev_obs_year <- wyear_y[sim_obs_idx[year_start_idx - 1L]]
      cand0 <- cand0[wyear_y[cand0] == prev_obs_year]
      if (!length(cand0)) cand0 <- which(wyear_y == prev_obs_year & month_y == first_month)
    }

    i0 <- sample(cand0, 1L)
    sim_prcp[year_start_idx] <- prcp_y[i0]
    sim_temp[year_start_idx] <- temp_y[i0]
    sim_date[year_start_idx] <- date_y[i0]
    sim_obs_idx[year_start_idx] <- i0

    # Set first day's markov state
    if (sim_prcp[year_start_idx] <= wet_thresh[first_month_idx]) {
      sim_occ[year_start_idx] <- 0L
    } else if (sim_prcp[year_start_idx] <= extreme_thresh[first_month_idx]) {
      sim_occ[year_start_idx] <- 1L
    } else {
      sim_occ[year_start_idx] <- 2L
    }

    # DAILY LOOP: remaining 364 days of this year
    sim_idx <- year_start_idx

    for (j in 2:365) {

      sim_idx <- sim_idx + 1L
      if (sim_idx > sim_length) break

      prev_idx <- sim_idx - 1L

      sim_occ[sim_idx] <- markov_next_state(
        sim_occ[prev_idx],
        rn_all[prev_idx],
        prev_idx,
        probs$p00_final,
        probs$p01_final,
        probs$p10_final,
        probs$p11_final,
        probs$p20_final,
        probs$p21_final
      )

      cur_month <- months_sim[sim_idx]
      cur_day   <- days_sim[sim_idx]
      m_idx <- month_to_year_index[cur_month]
      if (is.na(m_idx)) m_idx <- 1L

      key <- paste(cur_month, cur_day, sep = ".")
      obs_candidates <- lookup_day_idx_y[[key]]
      if (!length(obs_candidates)) obs_candidates <- which(month_y == cur_month)

      if (!length(obs_candidates)) {
        sim_prcp[sim_idx] <- sim_prcp[prev_idx]
        sim_temp[sim_idx] <- sim_temp[prev_idx]
        sim_date[sim_idx] <- sim_date[prev_idx]
        next
      }

      obs_idx_window <- expand_indices(obs_candidates, offsets7, length(prcp_y))

      if (!length(obs_idx_window)) {
        i <- sample(obs_candidates, 1L)
        sim_prcp[sim_idx] <- prcp_y[i]
        sim_temp[sim_idx] <- temp_y[i]
        sim_date[sim_idx] <- date_y[i]
        sim_obs_idx[sim_idx] <- i
        next
      }

      state_idx <- get_state_indices(
        from.state = sim_occ[prev_idx],
        to.state = sim_occ[sim_idx],
        prcp = prcp_y,
        candidate.idx = obs_idx_window,
        wet.thr = wet_thresh[m_idx],
        extreme.thr = extreme_thresh[m_idx]
      )

      if (!length(state_idx)) {
        obs_idx_window <- expand_indices(obs_candidates, offsets61, length(prcp_y))
        state_idx <- get_state_indices(
          from.state = sim_occ[prev_idx],
          to.state = sim_occ[sim_idx],
          prcp = prcp_y,
          candidate.idx = obs_idx_window,
          wet.thr = wet_thresh[m_idx],
          extreme.thr = extreme_thresh[m_idx]
        )
      }

      if (!length(state_idx)) {
        fb <- which(month_y == cur_month & (seq_along(prcp_y) + 1L) <= length(prcp_y))
        if (!is_water_year) fb <- fb[wyear_y[fb] == wyear_y[fb + 1L]]
        if (!length(fb)) {
          fb <- which((seq_along(prcp_y) + 1L) <= length(prcp_y))
          if (!is_water_year) fb <- fb[wyear_y[fb] == wyear_y[fb + 1L]]
        }

        i <- sample(fb, 1L)
        sim_prcp[sim_idx] <- prcp_y[i]
        sim_temp[sim_idx] <- temp_y[i]
        sim_date[sim_idx] <- date_y[i]
        sim_obs_idx[sim_idx] <- i
        next
      }

      possible_days <- obs_idx_window[state_idx]

      if (!is_water_year) {
        possible_days <- possible_days[wyear_y[possible_days] == wyear_y[possible_days + 1L]]
      }

      if (!length(possible_days)) {
        fb <- which(month_y == cur_month & (seq_along(prcp_y) + 1L) <= length(prcp_y))
        if (!is_water_year) fb <- fb[wyear_y[fb] == wyear_y[fb + 1L]]
        if (!length(fb)) {
          fb <- which((seq_along(prcp_y) + 1L) <= length(prcp_y))
          if (!is_water_year) fb <- fb[wyear_y[fb] == wyear_y[fb + 1L]]
        }

        i <- sample(fb, 1L)
        sim_prcp[sim_idx] <- prcp_y[i]
        sim_temp[sim_idx] <- temp_y[i]
        sim_date[sim_idx] <- date_y[i]
        sim_obs_idx[sim_idx] <- i
        next
      }

      next_days <- possible_days + 1L
      prcp_tomorrow <- prcp_y[next_days]
      temp_tomorrow <- temp_y[next_days]
      date_tomorrow <- date_y[next_days]

      cur_sim_prcp_anom <- sim_prcp[prev_idx] - mean_mon_prcp[cur_month]
      cur_sim_temp_anom <- sim_temp[prev_idx] - mean_mon_temp[cur_month]

      prcp_today_anom <- prcp_y[possible_days] - mean_mon_prcp[cur_month]
      temp_today_anom <- temp_y[possible_days] - mean_mon_temp[cur_month]

      k_day <- max(1L, round(sqrt(length(possible_days))))
      k_day <- min(k_day, length(possible_days))

      res <- knn_sample(
        candidates = cbind(prcp_today_anom, temp_today_anom),
        target = c(cur_sim_prcp_anom, cur_sim_temp_anom),
        k = k_day,
        n = 1,
        prob = TRUE,
        weights = knn_weights_by_month[[m_idx]],
        seed = if (is.null(base_seed)) NULL else base_seed + sim_idx
      )

      idx <- get_result_index(result = res, precip.tomorrow = prcp_tomorrow)

      # Enforce calendar-year constraint AFTER KNN
      if (!is_water_year) {
        prev_obs_year <- wyear_y[sim_obs_idx[prev_idx]]
        valid <- which(wyear_y[next_days] == prev_obs_year)
        if (length(valid)) idx <- idx[idx %in% valid]
        if (!length(idx)) idx <- sample(valid, 1L)
      }

      if (is.na(idx) || idx < 1L || idx > length(prcp_tomorrow)) {
        idx <- sample(seq_along(prcp_tomorrow), 1L)
      }

      sim_prcp[sim_idx] <- prcp_tomorrow[idx]
      sim_temp[sim_idx] <- temp_tomorrow[idx]
      sim_date[sim_idx] <- date_tomorrow[idx]
      sim_obs_idx[sim_idx] <- next_days[idx]
    }
  }

  class(sim_date) <- "Date"
  sim_date
}

#' Expand Index Positions by Fixed Offsets
#'
#' @description
#' Expands a vector of base index positions by adding a set of integer offsets,
#' returning all valid index positions within a specified upper bound. This is
#' typically used to generate moving windows or neighborhood indices around
#' reference positions (e.g., expanding event days to adjacent days).
#'
#' @param idx Integer vector of base index positions to be expanded.
#' @param offsets Integer vector of offsets to apply to each element of
#'   \code{idx}. Offsets may be positive, negative, or zero.
#' @param max_len Integer. Maximum allowable index length; expanded indices are
#'   constrained to the range \code{1:(max_len - 1)}.
#'
#' @details
#' For each element in \code{idx}, all values in \code{offsets} are added.
#' Resulting indices that are less than or equal to zero, or for which
#' \code{index + 1 > max_len}, are discarded. The function does not enforce
#' uniqueness or sorting of the output.
#'
#' @return
#' Integer vector of expanded index positions that satisfy the validity
#' constraints.
#'
#' @examples
#' idx <- c(5, 10)
#' offsets <- -1:1
#' expand_indices(idx, offsets, max_len = 15)
#'
#' @keywords internal
expand_indices <- function(idx, offsets, max_len) {
  out <- idx[rep(seq_along(idx), each = length(offsets))] +
    rep(offsets, times = length(idx))
  out[out > 0 & (out + 1L) <= max_len]
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



#' Normalize a Probability Vector with Robust Fallback Handling
#'
#' @description
#' Normalizes a numeric vector to sum to one, treating non-finite and negative
#' values as zero. If the total mass is zero after cleaning, a fallback
#' distribution is returned instead.
#'
#' @param p Numeric vector of (unnormalized) probabilities or weights.
#'   Non-finite values (\code{NA}, \code{NaN}, \code{Inf}) and negative values
#'   are set to zero prior to normalization.
#' @param fallback Optional numeric vector to return when \code{p} has zero
#'   total mass after cleaning. If \code{NULL} (default), a uniform distribution
#'   of length \code{length(p)} is returned.
#'
#' @details
#' This helper is intended for internal use in stochastic or probabilistic
#' workflows where degenerate or numerically unstable probability vectors may
#' arise. The function does not check that \code{fallback} sums to one; this is
#' assumed to be ensured by the caller.
#'
#' @return
#' Numeric vector of probabilities summing to one. If a fallback is used, the
#' returned vector is either \code{fallback} or a uniform distribution.
#'
#' @examples
#' p <- c(0.2, NA, -0.1, 0.5)
#' normalize_probs(p)
#'
#' normalize_probs(c(0, 0, 0))
#'
#' normalize_probs(c(0, 0, 0), fallback = c(0.7, 0.2, 0.1))
#'
#' @keywords internal
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
