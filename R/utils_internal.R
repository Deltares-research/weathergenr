

# ==============================================================================
# UTILS: STATISTICS
# ==============================================================================

#' Compute sample skewness
#'
#' @description
#' Computes the (non-bias-corrected) third standardized central moment (skewness) for a numeric
#' vector, excluding \code{NA} values.
#'
#' @param x Numeric vector.
#'
#' @details
#' Returns \code{NA} if fewer than 3 non-missing values are available or if the standard
#' deviation is zero (which yields undefined standardization).
#'
#' @return Numeric scalar. Skewness estimate, or \code{NA}.
#'
#' @keywords internal
compute_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)

  m <- mean(x)
  s <- sd(x)
  sum((x - m)^3) / (n * s^3)
}


#' Compute sample excess kurtosis
#'
#' @description
#' Computes the (non-bias-corrected) fourth standardized central moment minus 3 (excess kurtosis)
#' for a numeric vector, excluding \code{NA} values.
#'
#' @param x Numeric vector.
#'
#' @details
#' Returns \code{NA} if fewer than 4 non-missing values are available or if the standard
#' deviation is zero.
#'
#' @return Numeric scalar. Excess kurtosis estimate, or \code{NA}.
#'
#' @keywords internal
compute_kurtosis <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(NA)

  m <- mean(x)
  s <- sd(x)
  sum((x - m)^4) / (n * s^4) - 3  # Excess kurtosis
}


# ==============================================================================
# UTILS: RUNS / SPELLS
# ==============================================================================

#' Compute run-lengths (spell lengths) for dry or wet states
#'
#' @description
#' Computes consecutive run lengths of dry or wet states relative to a threshold.
#'
#' @param x Numeric vector. Precipitation values.
#' @param threshold Numeric scalar. Wet-day threshold.
#' @param below Logical. If \code{TRUE}, defines spells where \code{x < threshold} (dry spells).
#'   If \code{FALSE}, defines spells where \code{x >= threshold} (wet spells).
#'
#' @details
#' Uses \code{rle()} on the logical state series to identify consecutive runs. If no spells of the
#' requested type occur, the function returns \code{numeric(0)}.
#'
#' @return Numeric vector of positive integers giving spell lengths. May be length zero.
#'
#' @keywords internal
compute_spell_lengths <- function(x, threshold, below = TRUE) {

  if (below) {
    state <- x < threshold
  } else {
    state <- x >= threshold
  }

  # Find runs of TRUE
  rle.result <- rle(state)
  spell.lengths <- rle.result$lengths[rle.result$values]

  if (length(spell.lengths) == 0) {
    return(numeric(0))
  }

  spell.lengths
}


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


#' Assign qualitative assessments for moment preservation
#'
#' @description
#' Classifies moment changes based on absolute percent change thresholds that differ by metric.
#'
#' @param moments.df Data.frame. Output from \code{\link{compute_moment_diagnostics}} containing
#'   at least \code{metric} and \code{pct.change}.
#'
#' @details
#' The function applies metric-specific thresholds:
#' \itemize{
#'   \item mean/variance: excellent < 5, good < 15, else poor
#'   \item cv: excellent < 10, good < 20, else poor
#'   \item others (sd, skewness, kurtosis): good < 15, acceptable < 30, else poor
#' }
#'
#' @return Character vector of length \code{nrow(moments.df)} with assessment labels.
#'
#' @keywords internal
assess_moment_changes <- function(moments.df) {
  assessments <- character(nrow(moments.df))

  for (i in seq_len(nrow(moments.df))) {
    metric <- moments.df$metric[i]
    pct.change <- abs(moments.df$pct.change[i])

    if (metric %in% c("mean", "variance")) {
      assessments[i] <- ifelse(pct.change < 5, "excellent",
                               ifelse(pct.change < 15, "good", "poor"))
    } else if (metric == "cv") {
      assessments[i] <- ifelse(pct.change < 10, "excellent",
                               ifelse(pct.change < 20, "good", "poor"))
    } else {
      assessments[i] <- ifelse(pct.change < 15, "good",
                               ifelse(pct.change < 30, "acceptable", "poor"))
    }
  }

  assessments
}

# ==============================================================================
# UTILS: MISC
# ==============================================================================

#' Test Whether an Object Is a Finite Integer Scalar
#'
#' @description
#' Checks whether an object is a numeric scalar representing a finite integer
#' value. This helper is intended for lightweight input validation where strict
#' integer typing is not required but integer-valued numerics are acceptable.
#'
#' @param x Object to test.
#'
#' @details
#' The function returns \code{TRUE} if and only if:
#' \itemize{
#'   \item \code{x} is numeric,
#'   \item \code{x} has length 1,
#'   \item \code{x} is finite (not \code{NA}, \code{NaN}, or \code{Inf}),
#'   \item \code{x} has no fractional component (\code{x \%\% 1 == 0}).
#' }
#' Logical, character, and integer vectors of length greater than one will
#' return \code{FALSE}.
#'
#' @return
#' Logical scalar indicating whether \code{x} is a finite integer-valued scalar.
#'
#' @keywords internal
.is_int_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && (x %% 1 == 0)
}

#' Safely Compute Correlation with Pairwise Completeness Check
#'
#' @description
#' Computes a correlation coefficient between two numeric vectors using only
#' finite paired observations. If the number of valid pairs is below a required
#' minimum, the correlation is not computed.
#'
#' @param x,y Numeric vectors of equal length. Only pairs where both values are
#'   finite are used.
#' @param min_pairs Integer >= 1. Minimum number of paired observations required
#'   to compute the correlation.
#' @param method Character. Correlation method passed to \code{stats::cor()}.
#'   One of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#'
#' @return
#' Named numeric vector with elements:
#' \itemize{
#'   \item \code{value}: correlation coefficient, or \code{NA_real_},
#'   \item \code{n}: number of finite paired observations.
#' }
#'
#' @keywords internal
.safe_cor <- function(x, y, min_pairs = 3L, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)

  if (n < min_pairs) {
    return(c(value = NA_real_, n = n))
  }

  c(
    value = stats::cor(x[ok], y[ok], method = method),
    n     = n
  )
}




# ==============================================================================
# INTERNAL LOGGING UTILITIES
# ==============================================================================

#' Internal logging wrapper (INFO level)
#'
#' @description
#' Package-wide internal logging helper. Routes messages to the \pkg{logger}
#' package if available; otherwise falls back to \code{message()}.
#'
#' Logging is silent unless \code{verbose = TRUE}.
#'
#' @param msg Character scalar. Log message (already formatted).
#' @param verbose Logical. If FALSE, no output is produced.
#' @param tag Optional character scalar. Component tag (e.g. "FILTERING", "QMAP").
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
.log_info <- function(msg, verbose = FALSE, tag = NULL) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  if (!is.null(tag) && nzchar(tag)) {
    prefix <- paste0("[", tag, "]")
  }

  # Resolve { } HERE, in caller environment
  rendered <- msg
  if (requireNamespace("glue", quietly = TRUE)) {
    rendered <- tryCatch(
      as.character(glue::glue(msg, .envir = parent.frame())),
      error = function(e) msg
    )
  }

  full_msg <- paste(prefix, rendered)

  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(full_msg, .skip_formatter = TRUE)
  } else {
    message(full_msg)
  }

  invisible(NULL)
}


