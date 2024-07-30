
#' Calculate markov chain probobabilities
#'
#' @param p00 placeholder
#' @param p01 placeholder
#' @param p02 placeholder
#' @param p10 placeholder
#' @param p11 placeholder
#' @param p12 placeholder
#' @param p20 placeholder
#' @param p21 placeholder
#' @param p22 placeholder
#'
#' @keywords internal
#' @return a vector of markov probabilities
#' @export
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

  		PI <- solve(P[c(1,2,4),],B[c(1,2,4)])

  		return(PI)
	}
