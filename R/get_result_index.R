#' Safely Get or Sample an Index from a Vector
#'
#' Returns a valid index from `PRCP_TOMORROW` based on the provided `RESULT`.
#' If `RESULT` is missing, out of bounds, or invalid, a random index is sampled instead.
#' If `PRCP_TOMORROW` is empty, `NA_integer_` is returned.
#'
#' This function is useful in stochastic weather generation workflows where fallback
#' behavior is needed when a computed index is not valid.
#'
#' @param RESULT Integer. Proposed index for retrieving an element from `PRCP_TOMORROW`.
#' @param PRCP_TOMORROW Numeric vector. A set of candidate precipitation values (e.g., for the next day).
#'
#' @return Integer. A valid index from `PRCP_TOMORROW`, or `NA_integer_` if `PRCP_TOMORROW` is empty.
#'
#' @examples
#' PRCP_TOMORROW <- c(0, 5, 10, 20)
#'
#' # Valid index
#' get_result_index(2, PRCP_TOMORROW)  # returns 2
#'
#' # Invalid index: falls back to sampling
#' set.seed(1)
#' get_result_index(10, PRCP_TOMORROW)  # randomly returns a valid index
#'
#' # NA index: falls back to sampling
#' get_result_index(NA, PRCP_TOMORROW)
#'
#' # Empty vector: returns NA
#' get_result_index(1, numeric(0))
#'
#' @export
get_result_index <- function(RESULT, PRCP_TOMORROW) {
  if (is.na(RESULT) || length(RESULT) == 0 || RESULT < 1 || RESULT > length(PRCP_TOMORROW)) {
    if (length(PRCP_TOMORROW) > 0) sample(seq_along(PRCP_TOMORROW), 1)
    else NA_integer_
  } else {
    RESULT
  }
}
