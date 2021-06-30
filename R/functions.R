
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
