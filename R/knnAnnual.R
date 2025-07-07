
#' KNN-based Annual Weather Selection
#'
#' @param sim_annual_prcp Numeric vector of simulated annual precipitation
#' @param ANNUAL_PRCP Numeric vector of observed annual precipitation
#' @param WATER_YEAR_A Vector of water year indices or identifiers
#' @param kk Integer number of nearest neighbors to consider
#' @param k1 Integer, seed modifier (often a simulation run identifier)
#' @param y Integer, year index or current iteration
#' @param y_sample_size Integer, number of years to sample
#' @param seed Optional integer random seed for reproducibility
#'
#' @return
#'
#' Example: Select 5 nearest years for a simulated year
#' set.seed(42)
#' observed_prcp <- c(800, 850, 900, 780, 870, 820, 830, 890, 810, 860)
#' years <- 2001:2010
#' sim_year_prcp <- 845
#' sampled_years <- knnAnnual(
#'   sim_annual_prcp = sim_year_prcp,
#'   ANNUAL_PRCP = observed_prcp,
#'   WATER_YEAR_A = years,
#'   kk = 4,
#'   k1 = 2,
#'   y = 1,
#'   y_sample_size = 6,
#'   seed = 101
#' )
#' sampled_years  # Returns a vector of sampled years, e.g. c(850, 820, ...)
#'
#'
#' @export
knnAnnual <- function(sim_annual_prcp, ANNUAL_PRCP, WATER_YEAR_A,
                      kk, k1, y, y_sample_size = 50, seed = NULL) {

  # Handle random seed for reproducibility without global side-effects
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed + k1 * y)
  }

  # Calculate Euclidean distances
  distances <- sqrt((sim_annual_prcp - ANNUAL_PRCP)^2)

  # Identify k nearest neighbors
  k_distances <- order(distances)[1:kk]

  # Probability weighting for sampling
  if (kk == 1) {
    selection <- rep(k_distances, y_sample_size)
  } else {
    probs <- (1 / 1:kk) / sum(1 / 1:kk)
    selection <- sample(k_distances, y_sample_size, replace = TRUE, prob = probs)
  }

  # Return corresponding water year indices
  return(WATER_YEAR_A[selection])
}

