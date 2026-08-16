

# ==============================================================================
# UTILS: STATISTICS
# ==============================================================================

#' Compute sample skewness
#'
#' @description
#' Computes sample skewness as \eqn{\sum (x-\bar{x})^3 / (n s^3)} for a numeric
#' vector, excluding \code{NA} values.
#'
#' @param x Numeric vector.
#'
#' @details
#' Returns \code{NA} if fewer than 3 non-missing values are available or if the standard
#' deviation is zero (which yields undefined standardization). This estimator
#' uses \code{stats::sd()} (denominator \eqn{n-1}) in the standardization term
#' and follows the same convention used by \code{e1071::skewness(type = 2)}.
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
#' Computes sample excess kurtosis as
#' \eqn{\sum (x-\bar{x})^4 / (n s^4) - 3} for a numeric vector, excluding
#' \code{NA} values.
#'
#' @param x Numeric vector.
#'
#' @details
#' Returns \code{NA} if fewer than 4 non-missing values are available or if the standard
#' deviation is zero. As in \code{compute_skewness()}, the standardization uses
#' \code{stats::sd()} (denominator \eqn{n-1}).
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

#' Seed the RNG reproducibly, pinning the generator as well as the seed
#'
#' @description
#' A bare `set.seed(n)` keeps whatever generator is currently active, so the same
#' integer means different things in different contexts. A PSOCK worker prepared
#' by `parallel::clusterSetRNGStream()` runs L'Ecuyer-CMRG, while the master runs
#' Mersenne-Twister; `set.seed(12345)` in each then produces entirely different
#' streams. That is what made `generate_weather(parallel = TRUE)` disagree with
#' `parallel = FALSE` for one `seed`.
#'
#' Pinning the kind makes a seed mean one stream everywhere. It is a no-op in the
#' master, whose default generator is already Mersenne-Twister.
#'
#' Callers that save `.Random.seed` and restore it on exit also restore the
#' generator, because `.Random.seed[1]` encodes the kind -- so this does not leak
#' a generator change back to the caller.
#'
#' @param seed Integer scalar seed.
#' @return Invisibly `NULL`, called for its side effect.
#' @keywords internal
#' @noRd
.set_seed_fixed_kind <- function(seed) {
  # normal.kind and sample.kind are left NULL, i.e. unchanged.
  set.seed(seed, kind = "Mersenne-Twister")
  invisible(NULL)
}

#' Offset a seed without overflowing to NA
#'
#' @description
#' Derived seeds are built by adding a small offset to a base seed drawn from
#' `sample.int(.Machine$integer.max, 1L)`. Integer addition overflows to `NA`
#' when the base sits within the offset of the maximum, and `set.seed(NA)` is an
#' error -- a hard failure deep into a long run, reproducible only under the one
#' seed that triggers it.
#'
#' The arithmetic is done in double precision and wrapped, so every value that
#' did not overflow is returned unchanged. That matters: this is used on paths
#' the recorded numeric baseline depends on, and a fix that shifted ordinary
#' seeds would move output for no reason.
#'
#' Call a function, omitting arguments that are NULL
#'
#' @description
#' Forwards only the named arguments that are non-`NULL`, so an absent one falls
#' through to the callee's own default instead of overriding it.
#'
#' @details
#' `run_weather_generator()` builds its call from a user `config` list and
#' documents that "an entry that is absent (`NULL`) falls back to the receiving
#' function's default". Passing `config$x` directly does not do that: an absent
#' entry is an explicit `NULL`, which replaces the default and then fails the
#' callee's validation. A minimal config died on `dry_spell_factor must have
#' length 12` rather than defaulting to `rep(1, 12)`.
#'
#' Dropping the `NULL`s before the call is what makes the documented contract
#' true. Note this means a config cannot explicitly request `NULL` for an
#' argument whose default is non-`NULL` -- which is what "absent falls back to
#' the default" means, and none of the forwarded arguments treat `NULL` as a
#' meaningful value distinct from absence.
#'
#' @param .f Function to call.
#' @param ... Named arguments; those that are `NULL` are dropped.
#' @return The value of `.f` called with the retained arguments.
#' @keywords internal
#' @noRd
.call_dropping_null <- function(.f, ...) {
  args <- list(...)
  do.call(.f, args[!vapply(args, is.null, logical(1))])
}


#' @param seed Numeric scalar base seed.
#' @param offset Numeric scalar or vector added to `seed`.
#' @return Integer vector of seeds, wrapped into `[1, .Machine$integer.max]`.
#' @keywords internal
#' @noRd
.seed_offset <- function(seed, offset) {
  m <- .Machine$integer.max
  x <- as.numeric(seed) + as.numeric(offset)

  # Wrap into 1..m rather than 0..m-1 so a wrapped value is still a usable seed.
  as.integer(((x - 1) %% m) + 1)
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


#' Collapse a value that rounds to zero onto positive zero
#'
#' @description
#' `sprintf("%+.3f", -1e-16)` prints `-0.000`. The sign there is noise: a
#' residual that rounds to zero at the displayed precision may land on either
#' side of zero depending on the platform's floating point, so the printed sign
#' is neither stable nor meaningful. Snapping such values to `0` makes the
#' output identical across platforms and removes a nonsensical `-0%`.
#'
#' @param x Numeric vector.
#' @param digits Integer. Displayed precision the value is judged against.
#' @return `x` with values rounding to zero replaced by `0`.
#' @keywords internal
#' @noRd
.zap_signed_zero <- function(x, digits = 0) {
  ifelse(is.na(x) | round(x, digits) != 0, x, 0)
}

.format_pct <- function(x, digits = 0) {
  x <- .zap_signed_zero(x, digits)
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
#' - Lines are formatted as \code{HH:MM:SS - tag - message}, matching the
#'   downstream `blueearth_cst` console syntax
#'
#' @param msg Character scalar. Log message template with \code{{variable}} syntax.
#' @param level Character scalar. One of "info", "warn", "error".
#' @param verbose Logical. If FALSE, suppress output.
#' @param tag Optional character scalar. Component tag (e.g. "WARM", "KNN"),
#'   emitted lower-cased.
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

  # Coerce safely to a scalar character string
  if (length(msg) != 1L) msg <- paste(msg, collapse = " ")
  msg <- as.character(msg)

  # ---------------------------------------------------------------------------
  # Simple brace interpolation using base R (caller env)
  # ---------------------------------------------------------------------------
  rendered <- msg
  env <- parent.frame()

  matches <- gregexpr("\\{[^}]+\\}", msg)
  if (matches[[1]][1] != -1L) {
    exprs <- regmatches(msg, matches)[[1]]

    for (expr in exprs) {
      var_expr <- substr(expr, 2L, nchar(expr) - 1L)

      value <- tryCatch(
        eval(parse(text = var_expr), envir = env),
        error = function(e) expr
      )

      # Ensure replacement is a scalar string
      if (length(value) != 1L) value <- paste(value, collapse = ", ")
      value <- as.character(value)

      # Replace all occurrences of the exact token (fixed)
      rendered <- gsub(expr, value, rendered, fixed = TRUE)
    }
  }

  timestamp <- format(Sys.time(), "%H:%M:%S")
  if (!is.null(tag) && nzchar(tag)) {
    prefix <- paste0(timestamp, " - ", tolower(tag), " - ")
  } else {
    prefix <- paste0(timestamp, " - ")
  }

  rendered <- paste0(prefix, rendered)

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
