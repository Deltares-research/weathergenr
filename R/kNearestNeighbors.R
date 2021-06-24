

#' K nearest neighbor sampling
#'
#' This function applies a simple K nearest neighbors sampling scheme to find
#' the historical nearest neighbors of the current variable value and to resample
#' from k-nearest values. Adopted from Rajagopalan and Lall (1999)
#'
#' @param x A numeric vector of time-series of weather variables xt, t=1,..,n
#' @param k A numeric value to set the number of neighbors to resample from
#' @param s A numeric value for calculating the euclidean distance vector
#'
#' @return numeric, resampled value from the time-series vector x
#' @export
#'
kNearestNeighbors  <- function(x, k, s, seed = NULL) {

        dis <- sqrt((s - x)^2)
        s.ind <- order(dis)[1:k]
        s.wgh <- sapply(1:k, function(y) (1/y)/(sum(1/(1:k))))

        if(!is.null(seed)) {set.seed(seed)}
        val <- sample(s.ind, size = 1, replace = TRUE, prob = s.wgh)

        return(val)
}
