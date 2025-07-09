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
  if (prev_state == 0) {
    pp1 <- p00[idx]; pp2 <- p00[idx] + p01[idx]
  } else if (prev_state == 1) {
    pp1 <- p10[idx]; pp2 <- p10[idx] + p11[idx]
  } else {
    pp1 <- p20[idx]; pp2 <- p20[idx] + p21[idx]
  }
  if (is.na(pp1)) pp1 <- 0
  if (is.na(pp2)) pp2 <- 0
  if (rn < pp1) 0 else if (rn < pp2) 1 else 2
}
