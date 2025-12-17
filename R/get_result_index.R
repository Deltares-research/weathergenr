#' Safely Get or Sample an Index from a Vector
#'
#' Returns a valid index from `precip.tomorrow` based on the provided `result`.
#' If `result` is missing, out of bounds, or invalid, a random index is sampled instead.
#' If `precip.tomorrow` is empty, `NA_integer_` is returned.
#'
#' This function is useful in stochastic weather generation workflows where fallback
#' behavior is needed when a computed index is not valid.
#'
#' @param result Integer. Proposed index for retrieving an element from `precip.tomorrow`.
#' @param precip.tomorrow Numeric vector. A set of candidate precipitation values (e.g., for the next day).
#'
#' @return Integer. A valid index from `precip.tomorrow`, or `NA_integer_` if `precip.tomorrow` is empty.
#'
#' @examples
#' precip.tomorrow <- c(0, 5, 10, 20)
#'
#' # Valid index
#' get_result_index(2, precip.tomorrow) # returns 2
#'
#' # Invalid index: falls back to sampling
#' set.seed(1)
#' get_result_index(10, precip.tomorrow) # randomly returns a valid index
#'
#' # NA index: falls back to sampling
#' get_result_index(NA, precip.tomorrow)
#'
#' # Empty vector: returns NA
#' get_result_index(1, numeric(0))
#'
#' @export
get_result_index <- function(result, precip.tomorrow) {
  if (is.na(result) || length(result) == 0 || result < 1 || result > length(precip.tomorrow)) {
    if (length(precip.tomorrow) > 0) {
      sample(seq_along(precip.tomorrow), 1)
    } else {
      NA_integer_
    }
  } else {
    result
  }
}
