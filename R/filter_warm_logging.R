# ==============================================================================
# Logging Functions for filter_warm_simulations() - ENHANCED VERSION
# ==============================================================================

#' Log initial setup information
#'
#' @description
#' Displays general information at the start of filtering.
#'
#' @param n_obs Number of observations in observed series
#' @param n_sim Number of years in simulated series
#' @param n_realizations Number of realizations
#' @param sample_num Target number to sample
#' @param relax_priority Relaxation priority vector
#'
#' @keywords internal
#' @export
log_filtering_setup <- function(n_obs, n_sim, n_realizations, sample_num, relax_priority) {
  log_info(strrep("=", 75))
  log_info("FILTERING SETUP")
  log_info(strrep("=", 75))
  log_info(sprintf("Observed series: %d years", n_obs))
  log_info(sprintf("Simulated series: %d years x %d realizations", n_sim, n_realizations))
  log_info(sprintf("Target: Sample %d realizations from pool", sample_num))
  log_info(sprintf("Relaxation priority: %s", paste(relax_priority, collapse = " > ")))
  log_info(sprintf("  (Filters relax left-to-right: %s relaxes FIRST, %s relaxes LAST)",
                   relax_priority[1], relax_priority[length(relax_priority)]))
  log_info(strrep("=", 75))
  log_info("")
}

#' Log major step progress
#'
#' @description
#' Displays progress message for major computational steps.
#'
#' @param step_name Name of the step
#' @param details Optional details string
#'
#' @keywords internal
#' @export
log_step <- function(step_name, details = NULL) {
  msg <- if (!is.null(details)) {
    sprintf(" %s - %s", step_name, details)
  } else {
    sprintf(" %s", step_name)
  }
  log_info(msg)
}

#' Log filter iteration details with table format
#'
#' @description
#' Prints iteration diagnostics in table format for all iterations.
#'
#' @param iter Iteration number
#' @param passes List of pass vectors
#' @param pool Vector of pool indices
#' @param n_total Total number of realizations
#' @param target Target pool size
#' @param bounds Bounds environment or list
#' @param tail_metrics Tail metrics list
#' @param wavelet_active Logical
#' @param wavelet_pars Wavelet parameters list
#' @param note Optional note string
#'
#' @keywords internal
#' @export
log_filter_iteration <- function(iter, passes, pool, n_total, target, bounds,
                                 tail_metrics, wavelet_active, wavelet_pars,
                                 note = NULL) {

  active <- names(passes)[!vapply(passes, is.null, logical(1))]
  pool_size <- length(pool)
  pool_pct <- if (n_total > 0) 100 * pool_size / n_total else 0

  # Determine status
  if (pool_size >= target) {
    status <- "TARGET REACHED"
    status_icon <- "[OK]"
  } else {
    status <- "Need more candidates"
    status_icon <- "[>>]"
  }

  # Iteration header
  if (iter == 0L) {
    log_info(strrep("-", 75))
    log_info(sprintf("ITERATION %d - Initial Evaluation", iter))
  } else {
    log_info("")
    log_info(strrep("-", 75))
    log_info(sprintf("ITERATION %d - %s", iter, note))
  }
  log_info(strrep("-", 75))

  # Show filter table
  filter_order <- c("mean", "sd", "tail_low", "tail_high", "wavelet")
  show_filters <- intersect(filter_order, active)

  if (length(show_filters) > 0) {
    # Table header
    log_info(sprintf("%-12s %10s %8s  %-30s", "Filter", "Passed", "Rate", "Criteria"))
    log_info(strrep("-", 75))

    # Table rows
    for (nm in show_filters) {
      if (!is.null(passes[[nm]])) {
        n_pass <- sum(passes[[nm]])
        rate <- sprintf("%.1f%%", 100 * mean(passes[[nm]]))
        crit <- criteria_string_compact(nm, bounds, tail_metrics, wavelet_active, wavelet_pars)

        log_info(sprintf("%-12s %10d %8s  %-30s", nm, n_pass, rate, crit))
      }
    }
    log_info(strrep("-", 75))
  }

  # Status line
  log_info(sprintf("%s Pool: %d / %d (%.1f%%) | Need: %d | Status: %s",
                   status_icon, pool_size, n_total, pool_pct, target, status))

  if (pool_size >= target) {
    log_info(strrep("=", 75))
  }

  invisible(NULL)
}

#' Log final summary
#'
#' @description
#' Displays final filtering results.
#'
#' @param pool_size Final pool size
#' @param n_total Total realizations
#' @param n_sampled Number sampled
#' @param relaxation_level Relaxation level reached
#'
#' @keywords internal
#' @export
log_final_summary <- function(pool_size, n_total, n_sampled, relaxation_level) {
  log_info("")
  log_info(strrep("=", 75))
  log_info("FILTERING COMPLETE")
  log_info(strrep("=", 75))
  log_info(sprintf("Final pool: %d / %d realizations (%.1f%%)",
                   pool_size, n_total, 100 * pool_size / n_total))
  log_info(sprintf("Sampled: %d realizations", n_sampled))
  log_info(sprintf("Relaxation level: %s", relaxation_level))
  log_info(strrep("=", 75))
}

#' Compact criteria string for a filter
#'
#' @description
#' Creates a compact human-readable string describing current filter criteria.
#'
#' @param filter_name Character filter name
#' @param bounds Bounds environment or list
#' @param tail_metrics Tail metrics list
#' @param wavelet_active Logical
#' @param wavelet_pars Wavelet parameters
#'
#' @return Character string (compact)
#'
#' @keywords internal
#' @export
criteria_string_compact <- function(filter_name, bounds, tail_metrics,
                                    wavelet_active, wavelet_pars) {

  if (filter_name == "mean") {
    return(sprintf("tol = %.4f", bounds$mean))
  }

  if (filter_name == "sd") {
    return(sprintf("tol = %.4f", bounds$sd))
  }

  if (filter_name == "tail_low") {
    return(sprintf("p=%.2f, log.tol=%.4f", bounds$tail.low.p, bounds$tail.tol.log))
  }

  if (filter_name == "tail_high") {
    return(sprintf("p=%.2f, log.tol=%.4f", bounds$tail.high.p, bounds$tail.tol.log))
  }

  if (filter_name == "wavelet") {
    if (!isTRUE(wavelet_active)) return("inactive")
    return(sprintf("sig.frac >= %.2f", bounds$sig.frac))
  }

  "NA"
}

#' Format numeric with specified digits
#'
#' @param x Numeric value
#' @param digits Number of decimal places
#'
#' @return Character string
#'
#' @keywords internal
#' @export
fmt_num <- function(x, digits = 4L) {
  if (is.null(x) || length(x) == 0L || any(!is.finite(x))) return("NA")
  formatC(x, format = "f", digits = digits)
}

#' Format percentage
#'
#' @param x Numeric proportion (0 to 1)
#' @param digits Number of decimal places
#'
#' @return Character string with percent sign
#'
#' @keywords internal
#' @export
fmt_pct <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

#' Pad string to the right
#'
#' @param x Character string
#' @param width Target width
#'
#' @return Padded character string
#'
#' @keywords internal
#' @export
pad_right <- function(x, width) {
  x <- as.character(x)
  n <- nchar(x, type = "width")
  ifelse(n >= width, substr(x, 1, width), paste0(x, strrep(" ", width - n)))
}

#' Log info message with FILTERING prefix
#'
#' @description
#' Wrapper for logger::log_info with fallback to message().
#' Automatically adds [FILTERING] prefix to all messages.
#'
#' @param txt Character message
#'
#' @keywords internal
#' @export
log_info <- function(txt) {
  msg <- sprintf("[FILTERING] %s", txt)
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("{msg}")
  } else {
    message(msg)
  }
  invisible(NULL)
}
