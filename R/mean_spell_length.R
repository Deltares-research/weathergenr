#' @title Mean Spell Length Above or Below a Threshold
#'
#' @description
#' Computes the mean length of contiguous runs ("spells") in a numeric
#' time series after threshold-based classification.
#'
#' The input vector is first converted to a binary occurrence series using
#' the supplied threshold. Consecutive days in the same state are grouped
#' into spells, and the mean spell length is calculated for either the
#' below-threshold or above-threshold state.
#'
#' This function is typically used to diagnose wet or dry spell persistence
#' in daily precipitation or other hydro-meteorological time series and is
#' suitable for validating stochastic weather generators and Markov-chain
#' occurrence models.
#'
#' @param x Numeric vector. Daily values of a weather variable
#'   (for example precipitation or temperature). Missing values are not
#'   explicitly handled and should be removed beforehand if present.
#' @param threshold Numeric scalar. Threshold used to classify days into
#'   two states. Values less than or equal to the threshold are considered
#'   "below-threshold"; values above the threshold are considered
#'   "above-threshold".
#' @param below Logical. If TRUE (default), the mean length of
#'   below-threshold spells (for example dry spells) is returned.
#'   If FALSE, the mean length of above-threshold spells
#'   (for example wet spells) is returned.
#'
#' @return
#' A single numeric value giving the mean spell length (in days) for the
#' selected state. Returns:
#' \itemize{
#'   \item \code{NA_real_} if \code{x} has zero length,
#'   \item \code{0} if no spells of the selected type are present.
#' }
#'
#' @details
#' A spell is defined as a maximal sequence of consecutive days belonging
#' to the same threshold-defined state. Spell lengths are computed from
#' transitions in the binary occurrence series derived from \code{x}.
#'
#' The function treats values exactly equal to the threshold as
#' below-threshold. This convention should be kept consistent with other
#' occurrence or Markov-state definitions used in the analysis.
#'
#' @seealso
#' \code{\link{markov_next_state}}
#'
#' @export
#' @keywords internal
mean_spell_length <- function(x, threshold = 0, below = TRUE) {

  # Validate input
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }

  n <- length(x)
  if (n == 0L) {
    return(NA_real_)
  }

  # Convert to binary: 0 = dry, 1 = wet (or vice versa, depending on 'below')
  binary <- ifelse(x <= threshold, 0, 1)

  # Identify transitions
  change_points <- c(which(diff(binary) != 0), n)
  spell_lengths <- diff(c(0L, change_points))
  spell_types <- binary[change_points]

  # Select spell type
  selected_lengths <- spell_lengths[spell_types == if (below) 0 else 1]

  if (length(selected_lengths) == 0) {
    return(0)
  }

  # Calculate average length
  return(mean(selected_lengths))
}
