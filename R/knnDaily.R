
#' KNN daily sampling
#'
#' @param cur_sim_PRCP  placeholder
#' @param cur_sim_TEMP placeholder
#' @param PRCP_TODAY placeholder
#' @param TEMP_TODAY placeholder
#' @param DATE_TOMORROW placeholder
#' @param k placeholder
#' @param sd_monthly_PRCP placeholder
#' @param sd_monthly_TEMP placeholder
#' @param mean_monthly_PRCP placeholder
#' @param mean_monthly_TEMP placeholder
#' @param k1 placeholder
#' @param count placeholder
#'
#' @return
#' @export
knnDaily <- function(cur_sim_PRCP, cur_sim_TEMP, PRCP_TODAY, TEMP_TODAY,
    DATE_TOMORROW, k, sd_monthly_PRCP, sd_monthly_TEMP,
    mean_monthly_PRCP, mean_monthly_TEMP, k1, count) {

	  w_PRCP <- 100/sd_monthly_PRCP
		w_TEMP <- 10/sd_monthly_TEMP
		var_order <- 1:length(PRCP_TODAY)

		distance <- sqrt(w_PRCP*((cur_sim_PRCP-mean_monthly_PRCP) - (PRCP_TODAY-mean_monthly_PRCP))^2 + w_TEMP*((cur_sim_TEMP-mean_monthly_TEMP) - (TEMP_TODAY-mean_monthly_TEMP))^2)

		K_Distances <-	var_order[order(distance)[1:k]]
		probs <- (1/(1:k)) / sum((1/1:k))

		set.seed(k1*count)
		selection <- sample(1:k, size=1, prob=probs)

		return(DATE_TOMORROW[K_Distances[selection]])

	}
