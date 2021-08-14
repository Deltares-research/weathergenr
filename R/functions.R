
#' A wrapper to surppress the function output messages
#'
#' See Hadley Wickham http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#'
#' @param x An R function to surpress the output messages
#'
#' @return
#' @export
#'
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


#' Skewness function (from package e1071)
#'
#' @param x a numeric vector containing the values whose skewness is to be computed.
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @param type an integer between 1 and 3 selecting one of the algorithms for computing skewness detailed below.
#'
#' @return
#' @export
skewness <- function (x, na.rm = FALSE, type = 3)
{
    if (any(ina <- is.na(x))) {
        if (na.rm)
            x <- x[!ina]
        else return(NA)
    }
    if (!(type %in% (1:3)))
        stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
    if (type == 2) {
        if (n < 3)
            stop("Need at least 3 complete observations.")
        y <- y * sqrt(n * (n - 1))/(n - 2)
    }
    else if (type == 3)
        y <- y * ((1 - 1/n))^(3/2)
    y
}
