#' Get Indices for a Specific Occurrence State Transition
#'
#' Returns the indices in a precipitation time series (`PRCP`) where a specified
#' state transition occurs, based on wet/dry/extreme day thresholds.
#'
#' The occurrence states are defined as:
#' - 0: Dry day (precipitation <= `wet_thr`)
#' - 1: Wet day (precipitation > `wet_thr` and <= `extreme_thr`)
#' - 2: Extreme wet day (precipitation > `extreme_thr`)
#'
#' @param cur_occ Integer. The current occurrence state (0 = dry, 1 = wet, 2 = extreme).
#' @param next_occ Integer. The next occurrence state (0 = dry, 1 = wet, 2 = extreme).
#' @param PRCP Numeric vector of precipitation values.
#' @param cur_day Integer vector of indices representing the current day in the time series.
#' @param wet_thr Numeric. Threshold separating dry and wet days.
#' @param extreme_thr Numeric. Threshold above which days are considered extreme.
#'
#' @return Integer vector of indices where the given transition from `cur_occ` to `next_occ` occurs.
#'
#' @examples
#' PRCP <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
#' cur_day <- 1:(length(PRCP) - 1)
#' wet_thr <- 1
#' extreme_thr <- 20
#'
#' # Find transitions from dry to wet
#' get_state_indices(0, 1, PRCP, cur_day, wet_thr, extreme_thr)
#'
#' # Find transitions from wet to extreme
#' get_state_indices(1, 2, PRCP, cur_day, wet_thr, extreme_thr)
#'
#' # Find transitions from extreme to dry
#' get_state_indices(2, 0, PRCP, cur_day, wet_thr, extreme_thr)
#'
#' @export
get_state_indices <- function(cur_occ, next_occ, PRCP, cur_day, wet_thr, extreme_thr) {
  if (cur_occ == 0 && next_occ == 0)
    which(PRCP[cur_day] <= wet_thr & PRCP[(cur_day + 1)] <= wet_thr)
  else if (cur_occ == 0 && next_occ == 1)
    which(PRCP[cur_day] <= wet_thr & PRCP[(cur_day + 1)] > wet_thr & PRCP[(cur_day + 1)] <= extreme_thr)
  else if (cur_occ == 0 && next_occ == 2)
    which(PRCP[cur_day] <= wet_thr & PRCP[(cur_day + 1)] > extreme_thr)
  else if (cur_occ == 1 && next_occ == 0)
    which(PRCP[cur_day] > wet_thr & PRCP[cur_day] <= extreme_thr & PRCP[(cur_day + 1)] <= wet_thr)
  else if (cur_occ == 1 && next_occ == 1)
    which(PRCP[cur_day] > wet_thr & PRCP[cur_day] <= extreme_thr & PRCP[(cur_day + 1)] > wet_thr & PRCP[(cur_day + 1)] <= extreme_thr)
  else if (cur_occ == 1 && next_occ == 2)
    which(PRCP[cur_day] > wet_thr & PRCP[cur_day] <= extreme_thr & PRCP[(cur_day + 1)] > extreme_thr)
  else if (cur_occ == 2 && next_occ == 0)
    which(PRCP[cur_day] > extreme_thr & PRCP[(cur_day + 1)] <= wet_thr)
  else if (cur_occ == 2 && next_occ == 1)
    which(PRCP[cur_day] > extreme_thr & PRCP[(cur_day + 1)] > wet_thr & PRCP[(cur_day + 1)] <= extreme_thr)
  else if (cur_occ == 2 && next_occ == 2)
    which(PRCP[cur_day] > extreme_thr & PRCP[(cur_day + 1)] > extreme_thr)
  else
    integer(0)
}
