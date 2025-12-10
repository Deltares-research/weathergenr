
#' KNN daily sampling
#'
#' @param cur_sim_PRCP  placeholder
#' @param cur_sim_TEMP placeholder
#' @param PRCP_TODAY placeholder
#' @param TEMP_TODAY placeholder
#' @param k placeholder
#' @param sd_monthly_PRCP placeholder
#' @param sd_monthly_TEMP placeholder
#' @param mean_monthly_PRCP placeholder
#' @param mean_monthly_TEMP placeholder
#' @param k1 placeholder
#' @param count placeholder
#' @param seed random seed value
#'
#' @return
#' PLACEHOLDER
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

  # If ne seed provided, generate a random number
  if(is.null(seed)) seed = sample.int(1e10,1)

  w_PRCP <- 100/sd_monthly_PRCP
	w_TEMP <- 10/sd_monthly_TEMP
	var_order <- 1:length(PRCP_TODAY)

	distance <- sqrt(w_PRCP*((cur_sim_PRCP-mean_monthly_PRCP) - (PRCP_TODAY-mean_monthly_PRCP))^2 +
	                 w_TEMP*((cur_sim_TEMP-mean_monthly_TEMP) - (TEMP_TODAY-mean_monthly_TEMP))^2)

	K_Distances <-	var_order[order(distance)[1:k]]

	set.seed(seed + k1*count)
	selection <- sample(1:k, 1, prob=(1/1:k)/sum(1/1:k))

	return(K_Distances[selection])

	}

