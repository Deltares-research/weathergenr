
#' Calculate Spell Length
#'
#' @param x  placeholder
#' @param threshold placeholder
#' @param below placeholder
#'
#' @return
#' @export
calculateSpellLength <- function (x, threshold = 0, below = TRUE)
{

    if (!is.vector(x) && !is.list(x))
        stop("'x' must be a vector of an atomic type")
    n <- length(x)

    if (n == 0L)
        return(structure(list(lengths = integer(), values = x)))

    #Translate vector to binary based on threshold
    x[x<=threshold] = 0
    x[x>threshold] = 1

    y <- x[-1L] != x[-n]
    i <- c(which(y | is.na(y)), n)

    if(isTRUE(below)) {
      diff(c(0L, i))[which(x[i] == 0)]
    } else {
      diff(c(0L, i))[which(x[i] == 1)]
    }
}

