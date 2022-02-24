
#' Calculate Spell Length
#'
#' @param x  placeholder
#' @param threshold placeholder
#' @param below placeholder
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

