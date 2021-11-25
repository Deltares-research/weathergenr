
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
  y_sample_size = 20)

  {

	  var_order <- 1:length(ANNUAL_PRCP)
		distance <- sqrt((sim_annual_prcp - ANNUAL_PRCP)^2)

		ordered_distances <- matrix(cbind(var_order,distance)[order(distance),],ncol=2)
		K_Distances <- matrix(ordered_distances[1:kk,],ncol=2)
		probs <- (1/row(K_Distances)[,1]) / sum(1/row(K_Distances)[,1])

		set.seed(k1*y)
		selection <- sample(row(K_Distances)[,1], y_sample_size, replace=TRUE, prob=probs)

		return(WATER_YEAR_A[K_Distances[selection,1]])
	}

