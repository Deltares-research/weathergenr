

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
  if (!is.finite(s) || s == 0) return(NA)
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
  if (!is.finite(s) || s == 0) return(NA)
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
    state <- x <= threshold
  } else {
    state <- x > threshold
  }

  rle_result <- rle(state)
  spell_lengths <- rle_result$lengths[rle_result$values]

  if (length(spell_lengths) == 0) return(numeric(0))
  spell_lengths
}





#' Assign qualitative assessments for moment preservation
#'
#' @description
#' Classifies moment changes based on absolute percent change thresholds that differ by metric.
#'
#' @param moments_df Data.frame. Output from \code{\link{compute_moment_diagnostics}} containing
#'   at least \code{metric} and \code{pct_change}.
#'
#' @details
#' The function applies metric-specific thresholds:
#' \itemize{
#'   \item mean/variance: excellent < 5, good < 15, else poor
#'   \item cv: excellent < 10, good < 20, else poor
#'   \item others (sd, skewness, kurtosis): good < 15, acceptable < 30, else poor
#' }
#'
#' @return Character vector of length \code{nrow(moments_df)} with assessment labels.
#'
#' @keywords internal
assess_moment_changes <- function(moments_df) {
  assessments <- character(nrow(moments_df))

  for (i in seq_len(nrow(moments_df))) {
    metric <- moments_df$metric[i]
    pct_change <- abs(moments_df$pct_change[i])

    if (metric %in% c("mean", "variance")) {
      assessments[i] <- ifelse(pct_change < 5, "excellent",
                               ifelse(pct_change < 15, "good", "poor"))
    } else if (metric == "cv") {
      assessments[i] <- ifelse(pct_change < 10, "excellent",
                               ifelse(pct_change < 20, "good", "poor"))
    } else {
      assessments[i] <- ifelse(pct_change < 15, "good",
                               ifelse(pct_change < 30, "acceptable", "poor"))
    }
  }

  assessments
}


#' Compute Area-Averaged Daily and Annual Climate Series
#'
#' @description
#' Computes area-averaged (mean across grid cells) daily climate values and
#' aggregates them to annual means by water year.
#'
#' @param obs_data Named list of data frames (one per grid cell).
#' @param wyear_idx Integer vector of row indices to extract from each data frame.
#' @param wyear Integer vector of water years corresponding to wyear_idx.
#' @param vars Character vector of variable names to average.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{daily}{Data frame of area-averaged daily values with columns for
#'     each variable plus wyear.}
#'   \item{annual}{Tibble of annual means with columns wyear plus each variable.}
#' }
#'
#' @export
compute_area_averages <- function(obs_data, wyear_idx, wyear, vars) {

  n_grids <- length(obs_data)
  n_days  <- length(wyear_idx)
  n_vars  <- length(vars)

  # ---------------------------------------------------------------------------
  # Daily area average
  # ---------------------------------------------------------------------------
  if (n_grids == 1L) {
    daily_avg <- obs_data[[1]][wyear_idx, vars, drop = FALSE]
  } else {
    daily_mat <- matrix(0, nrow = n_days, ncol = n_vars)
    colnames(daily_mat) <- vars

    for (i in seq_len(n_grids)) {
      grid_data <- as.matrix(obs_data[[i]][wyear_idx, vars, drop = FALSE])
      daily_mat <- daily_mat + grid_data
    }

    daily_avg <- as.data.frame(daily_mat / n_grids)
  }

  daily_avg$wyear <- wyear

  # ---------------------------------------------------------------------------
  # Annual area average (mean of daily values by water year)
  # ---------------------------------------------------------------------------
  annual_avg <- daily_avg |>
    dplyr::group_by(wyear) |>
    dplyr::summarize(dplyr::across(dplyr::all_of(vars), mean), .groups = "drop")

  list(
    daily  = daily_avg,
    annual = annual_avg
  )
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


.format_pct <- function(x, digits = 0) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%+.", digits, "f%%"), x))
}

.format_num <- function(x, digits = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}



# ==============================================================================
# INTERNAL LOGGING UTILITIES
# ==============================================================================
#' Package-wide internal logger
#'
#' @description
#' Unified internal logging helper for the package.
#'
#' Features:
#' - Single entry point for all logging
#' - Brace interpolation resolved in caller environment (base R, no glue)
#' - Supports log levels (info, warn, error)
#' - Silent unless verbose = TRUE
#' - Timestamps in ISO format
#'
#' @param msg Character scalar. Log message template with \code{{variable}} syntax.
#' @param level Character scalar. One of "info", "warn", "error".
#' @param verbose Logical. If FALSE, suppress output.
#' @param tag Optional character scalar. Component tag (e.g. "WARM", "KNN").
#'
#' @return Invisibly returns NULL.
#'
#' @keywords internal
.log <- function(msg,
                 level = c("info", "warn", "error"),
                 verbose = TRUE,
                 tag = NULL) {

  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  level <- match.arg(level)

  # ---------------------------------------------------------------------------
  # Simple brace interpolation using base R
  # ---------------------------------------------------------------------------
  rendered <- msg
  env <- parent.frame()

  # Find all {expression} patterns
  matches <- gregexpr("\\{[^}]+\\}", msg)
  if (matches[[1]][1] != -1L) {
    exprs <- regmatches(msg, matches)[[1]]
    for (expr in exprs) {
      var_expr <- substr(expr, 2L, nchar(expr) - 1L)
      value <- tryCatch(
        eval(parse(text = var_expr), envir = env),
        error = function(e) expr
      )
      rendered <- sub(expr, as.character(value), rendered, fixed = TRUE)
    }
  }

  # Build prefix: timestamp + optional tag
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

  if (!is.null(tag) && nzchar(tag)) {
    prefix <- paste0(timestamp, " [", tag, "] ")
  } else {
    prefix <- paste0(timestamp, " ")
  }

  rendered <- paste0(prefix, rendered)

  # ---------------------------------------------------------------------------
  # Emit via base R
  # ---------------------------------------------------------------------------
  switch(
    level,
    info  = message(rendered),
    warn  = warning(rendered, call. = FALSE),
    error = stop(rendered, call. = FALSE)
  )

  invisible(NULL)
}


#' Format Elapsed Time for Display
#'
#' @description
#' Computes elapsed time from a start time and formats it as a human-readable
#' string with appropriate units (seconds, minutes, or hours).
#'
#' @param start_time POSIXct object. The start time to measure from.
#'
#' @return Character string with formatted elapsed time.
#'
#' @keywords internal
format_elapsed <- function(start_time) {
  elapsed_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (elapsed_secs < 60) {
    sprintf("%.1f seconds", elapsed_secs)
  } else if (elapsed_secs < 3600) {
    mins <- floor(elapsed_secs / 60)
    secs <- round(elapsed_secs %% 60, 0)
    sprintf("%d min %d sec", mins, secs)
  } else {
    hrs <- floor(elapsed_secs / 3600)
    mins <- round((elapsed_secs %% 3600) / 60, 0)
    sprintf("%d hr %d min", hrs, mins)
  }
}


