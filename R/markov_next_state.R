#' Determine Next State in a First-Order Markov Chain
#'
#' Given the previous state and a random number `rn`, this function returns the next
#' state in a 3-state first-order Markov chain. The transition probabilities are
#' state- and index-specific.
#'
#' States are typically defined as:
#' - 0: Dry
#' - 1: Wet
#' - 2: Extreme
#'
#' The transition is determined using cumulative probabilities for the transitions
#' to state 0 and 1. Any remaining probability is assigned to state 2.
#'
#' @param prev_state Integer. The current (previous) state (0, 1, or 2).
#' @param rn Numeric. A random number in [0, 1] used to sample the next state.
#' @param idx Integer. Index for selecting transition probabilities (e.g., time step or spatial location).
#' @param p00 Numeric vector. Probability of transition from state 0 to 0.
#' @param p01 Numeric vector. Probability of transition from state 0 to 1.
#' @param p10 Numeric vector. Probability of transition from state 1 to 0.
#' @param p11 Numeric vector. Probability of transition from state 1 to 1.
#' @param p20 Numeric vector. Probability of transition from state 2 to 0.
#' @param p21 Numeric vector. Probability of transition from state 2 to 1.
#'
#' @return Integer. The next state (0, 1, or 2).
#'
#' @examples
#' # Suppose we're at time step 5 with current state = 1
#' set.seed(123)
#' rn <- runif(1)
#' markov_next_state(
#'   prev_state = 1,
#'   rn = rn,
#'   idx = 5,
#'   p00 = rep(0.7, 10),
#'   p01 = rep(0.2, 10),
#'   p10 = rep(0.3, 10),
#'   p11 = rep(0.4, 10),
#'   p20 = rep(0.1, 10),
#'   p21 = rep(0.3, 10)
#' )
#'
#' @export
markov_next_state <- function(prev_state, rn, idx, p00, p01, p10, p11, p20, p21) {

  # Clamp random number
  if (!is.finite(rn)) rn <- 0
  if (rn < 0) rn <- 0
  if (rn > 1) rn <- 1

  # Validate prev_state
  if (!(prev_state %in% c(0L, 1L, 2L))) {
    prev_state <- 0L
  }

  # Clamp idx to valid range
  n <- length(p00)
  if (idx < 1L) idx <- 1L
  if (idx > n)  idx <- n

  # Select row probabilities
  if (prev_state == 0L) {
    a <- p00[idx]
    b <- p01[idx]
  } else if (prev_state == 1L) {
    a <- p10[idx]
    b <- p11[idx]
  } else {
    a <- p20[idx]
    b <- p21[idx]
  }

  # NA -> 0
  if (!is.finite(a)) a <- 0
  if (!is.finite(b)) b <- 0

  # Clamp negatives
  if (a < 0) a <- 0
  if (b < 0) b <- 0

  # Enforce stochastic validity
  s <- a + b
  if (s > 1) {
    a <- a / s
    b <- b / s
  }

  # Draw next state
  if (rn < a) {
    0L
  } else if (rn < (a + b)) {
    1L
  } else {
    2L
  }
}
