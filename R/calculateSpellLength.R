
#' Calculate length of a wet or a dry spell
#'
#' @param x numeric vector of values
#' @param threshold numeric threshold to define the wet or a dry spell
#' @param below binary variable. TRUE considers values less than or equal, FALSE above the threshold.
#'
#' @return
#' @export
averageSpellLength <- function (x, threshold = 0, below = TRUE)
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
      res <- diff(c(0L, i))[which(x[i] == 0)]
    } else {
      res <- diff(c(0L, i))[which(x[i] == 1)]
    }

    spell_length = as.numeric(names(table(res)))
    spell_count = as.numeric(table(res))

    sum(spell_length * spell_count)/sum(spell_count)

}

