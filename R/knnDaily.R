#' @title K-Nearest Neighbors (KNN) Daily Weather Resampling
#'
#' @description
#' Selects the index of a historical day whose weather is most similar to a simulated day's weather,
#' using the K-nearest neighbors algorithm. Similarity is determined based on standardized precipitation
#' and temperature. The function allows for stochastic selection among the K closest matches, favoring
#' closer neighbors.
#'
#' @param cur_sim_PRCP Numeric scalar. Simulated precipitation value for the current day.
#' @param cur_sim_TEMP Numeric scalar. Simulated temperature value for the current day.
#' @param PRCP_TODAY Numeric vector. Historical daily precipitation values to compare against.
#' @param TEMP_TODAY Numeric vector. Historical daily temperature values to compare against.
#' @param k Integer. Number of nearest neighbors to consider (must be at least 1).
#' @param w_PRCP Numeric scalar. Weight for precipitation in distance calculation. Typically `100 / sd_monthly_PRCP`.
#' @param w_TEMP Numeric scalar. Weight for temperature in distance calculation. Typically `10 / sd_monthly_TEMP`.
#' @param seed Optional integer. If provided, ensures reproducibility of random selection among neighbors.
#'
#' @details
#' For each simulated day, the Euclidean distance between the simulated and historical days is computed
#' using weighted precipitation and temperature differences. Among the K nearest days, a single day is
#' stochastically selected with a probability inversely proportional to its rank (closer neighbors are
#' more likely).
#'
#' The weights `w_PRCP` and `w_TEMP` are typically set as the inverse of the standard deviations of precipitation
#' and temperature, scaled to balance their influence in the distance metric.
#'
#' If `k` exceeds the number of available historical days, it is automatically reduced.
#'
#' @return
#' Integer. Index of the selected historical day from `PRCP_TODAY`/`TEMP_TODAY`.
#' Returns `NA_integer_` if no valid selection is possible.
#'
#' @examples
#' set.seed(42)
#' knnDaily(
#'   cur_sim_PRCP = 5,
#'   cur_sim_TEMP = 22,
#'   PRCP_TODAY = runif(100, 0, 20),
#'   TEMP_TODAY = runif(100, 15, 30),
#'   k = 10,
#'   w_PRCP = 20,
#'   w_TEMP = 5,
#'   seed = 123
#' )
#'
#' @export
knnDaily <- function(
  cur_sim_PRCP = NULL,
  cur_sim_TEMP = NULL,
  PRCP_TODAY = NULL,
  TEMP_TODAY = NULL,
  k = NULL,
  w_PRCP = 100/sd_monthly_PRCP,
  w_TEMP = 10/sd_monthly_TEMP,
  seed = NULL)

{

  # # Defensive checks
  # if (length(PRCP_TODAY) != length(TEMP_TODAY)) {
  #   stop("PRCP_TODAY and TEMP_TODAY must have the same length.")
  # }
  # if (k > length(PRCP_TODAY)) {
  #   stop("k cannot be greater than the number of available days.")
  # }
  # if (anyNA(c(cur_sim_PRCP, cur_sim_TEMP, PRCP_TODAY, TEMP_TODAY,
  #             sd_monthly_PRCP, sd_monthly_TEMP, mean_monthly_PRCP, mean_monthly_TEMP))) {
  #   stop("Inputs contain NA values.")
  # }

  # Local random seed
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed)
  }

	# Never select more than available days
	k <- min(k, length(PRCP_TODAY))

	# Calculate distances
	distance <- sqrt(w_PRCP*(cur_sim_PRCP - PRCP_TODAY)^2 + w_TEMP*(cur_sim_TEMP - TEMP_TODAY)^2)

	var_order <- seq_along(PRCP_TODAY)
	k_distances <-	var_order[order(distance)[1:k]]
  if(k == 1) return(k_distances[1])

	# Defensive KNN selection
	if (length(k_distances) == 0 || is.null(k) || k == 0) {
	  return(NA_integer_)
	}

	# Probability weighting
	probs <- (1 / seq_len(k)) / sum(1 / seq_len(k))

	# Defensive: no NA/NaN, finite, not all zero
	# if (anyNA(probs) || any(!is.finite(probs)) || all(probs == 0)) {
	#   return(NA_integer_)
	# }

  selection <- sample(1:k, 1, prob = probs)

	# if (length(selection) == 0 || is.na(selection) || selection < 1 || selection > length(k_distances)) {
	#   return(NA_integer_)
	# }
	return(k_distances[selection])

}
