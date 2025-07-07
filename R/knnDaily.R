#' @title K-Nearest Neighbor (KNN) Daily Weather Resampling
#'
#' @description
#' Selects a historical day most similar to the current simulated weather using
#' K-nearest neighbor resampling based on standardized precipitation and temperature.
#'
#' @param cur_sim_PRCP Numeric value of current simulated precipitation.
#' @param cur_sim_TEMP Numeric value of current simulated temperature.
#' @param PRCP_TODAY Numeric vector of historical daily precipitation values.
#' @param TEMP_TODAY Numeric vector of historical daily temperature values.
#' @param k Number of nearest neighbors to consider.
#' @param sd_monthly_PRCP Monthly standard deviation of historical precipitation.
#' @param sd_monthly_TEMP Monthly standard deviation of historical temperature.
#' @param mean_monthly_PRCP Monthly mean of historical precipitation.
#' @param mean_monthly_TEMP Monthly mean of historical temperature.
#' @param k1 Integer index for simulation trace (used for seed offset).
#' @param count Integer time step (e.g., day index in simulation; used for seed offset).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An integer index corresponding to the selected historical day from `PRCP_TODAY`.
#'
#' @examples
#' set.seed(42)
#' knnDaily(
#'   cur_sim_PRCP = 5,
#'   cur_sim_TEMP = 22,
#'   PRCP_TODAY = runif(100, 0, 20),
#'   TEMP_TODAY = runif(100, 15, 30),
#'   k = 10,
#'   sd_monthly_PRCP = 5,
#'   sd_monthly_TEMP = 3,
#'   mean_monthly_PRCP = 10,
#'   mean_monthly_TEMP = 20,
#'   k1 = 1,
#'   count = 7,
#'   seed = 123)
#'
#' @export
knnDaily <- function(
  cur_sim_PRCP = NULL,
  cur_sim_TEMP = NULL,
  PRCP_TODAY = NULL,
  TEMP_TODAY = NULL,
  k = NULL,
  sd_monthly_PRCP = NULL,
  sd_monthly_TEMP = NULL,
  mean_monthly_PRCP = NULL,
  mean_monthly_TEMP = NULL,
  k1 = NULL,
  count  = NULL,
  seed = NULL)

{

  # Defensive checks
  if (length(PRCP_TODAY) != length(TEMP_TODAY)) {
    stop("PRCP_TODAY and TEMP_TODAY must have the same length.")
  }
  if (k > length(PRCP_TODAY)) {
    stop("k cannot be greater than the number of available days.")
  }
  if (anyNA(c(cur_sim_PRCP, cur_sim_TEMP, PRCP_TODAY, TEMP_TODAY,
              sd_monthly_PRCP, sd_monthly_TEMP, mean_monthly_PRCP, mean_monthly_TEMP))) {
    stop("Inputs contain NA values.")
  }

  # Local random seed
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed + k1 * count)
  }

  # Set weights
  w_PRCP <- 100/sd_monthly_PRCP
	w_TEMP <- 10/sd_monthly_TEMP
	var_order <- seq_along(PRCP_TODAY)

	# Never select more than available days
	k <- min(k, length(PRCP_TODAY))

	# Calculate distances
	distance <- sqrt(w_PRCP*((cur_sim_PRCP-mean_monthly_PRCP) - (PRCP_TODAY-mean_monthly_PRCP))^2 +
	                 w_TEMP*((cur_sim_TEMP-mean_monthly_TEMP) - (TEMP_TODAY-mean_monthly_TEMP))^2)

	k_distances <-	var_order[order(distance)[1:k]]

	# Defensive KNN selection
	if (length(k_distances) == 0 || is.null(k) || k == 0) {
	  return(NA_integer_)
	}

	# Probability weighting
	probs <- (1 / seq_len(k)) / sum(1 / seq_len(k))

	# Defensive: no NA/NaN, finite, not all zero
	if (anyNA(probs) || any(!is.finite(probs)) || all(probs == 0)) {
	  return(NA_integer_)
	}

	if (k == 1) {
	  selection <- 1
	} else {
	  selection <- sample(1:k, 1, prob = probs)
	}

	if (length(selection) == 0 || is.na(selection) || selection < 1 || selection > length(k_distances)) {
	  return(NA_integer_)
	}
	return(k_distances[selection])

}
