
#' This KNN function generates 100 years.....
#'
#' @param sim_annual_prcp placeholder
#' @param ANNUAL_PRCP placeholder
#' @param WATER_YEAR_A placeholder
#' @param kk placeholder
#' @param k1 placeholder
#' @param y placeholder
#' @param y_sample_size placeholder
#'
#' @return
#' @export
knnAnnual <- function(
  sim_annual_prcp = NULL,
  ANNUAL_PRCP = NULL,
  WATER_YEAR_A = NULL,
  kk = NULL,
	k1 = NULL,
  y = NULL,
  y_sample_size = 50)

  {

		distances <- sqrt((sim_annual_prcp - ANNUAL_PRCP)^2)

		k_distances <- order(distances)[1:kk]
		probs <- (1/1:kk) /sum(1/1:kk)

		set.seed(k1*y)
		selection <- sample(1:kk, y_sample_size, replace=TRUE, prob=probs)

		return(WATER_YEAR_A[k_distances[selection]])
	}

