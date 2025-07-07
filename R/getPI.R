
#' @title Compute Stationary Distribution of a 3-State Markov Chain
#'
#' @description
#' Solves for the stationary probabilities of a 3-state Markov chain
#' given the transition probabilities between all states.
#'
#' @param p00 Probability of transitioning from state 0 to 0
#' @param p01 Probability of transitioning from state 0 to 1
#' @param p02 Probability of transitioning from state 0 to 2
#' @param p10 Probability of transitioning from state 1 to 0
#' @param p11 Probability of transitioning from state 1 to 1
#' @param p12 Probability of transitioning from state 1 to 2
#' @param p20 Probability of transitioning from state 2 to 0
#' @param p21 Probability of transitioning from state 2 to 1
#' @param p22 Probability of transitioning from state 2 to 2
#'
#' @return A numeric vector of stationary probabilities.
#' @export
#' @keywords internal
getPI <- function(
  p00, p01, p02,
  p10, p11, p12,
  p20, p21, p22)

  {

	  pi0 <- 1
		pi1 <- 1
		pi2 <- 1
		P <- matrix(c(p00,p01,p02,pi0,p10,p11,p12,pi1,p20,p21,p22,pi2),byrow=FALSE,nrow=4)
		P[1,1] <- P[1,1] - 1
		P[2,2] <- P[2,2] - 1
		P[3,3] <- P[3,3] - 1
		B <- c(0,0,0,1)

		PI <- tryCatch({
		  solve(P[c(1,2,4),], B[c(1,2,4)])
		}, error = function(e) {
		  stop("Stationary distribution could not be solved: ", conditionMessage(e))
		})
		if (!(all(is.finite(PI)) && all(PI >= 0) && all(PI <= 1) && abs(sum(PI) - 1) < 1e-6)) {
		  stop("Stationary distribution is not valid (check transition matrix input).")
		}
		return(PI)

}

