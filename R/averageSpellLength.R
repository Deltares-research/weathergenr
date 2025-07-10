#' @title Calculate Average Spell Length
#' @description
#' Computes the average length of either wet or dry spells based on a numeric input vector and a threshold.
#'
#' @param x A numeric vector representing daily weather values (e.g., precipitation).
#' @param threshold A numeric threshold used to define spell classification.
#' @param below Logical. If TRUE, calculates dry spells (values less than threshold). If FALSE, calculates wet spells (values more than threshold).
#'
#' @return A single numeric value representing the average spell length.
#' @export
#' @keywords internal

averageSpellLength <- function(x, threshold = 0, below = TRUE) {
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
