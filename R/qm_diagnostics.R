# ==============================================================================
# VALIDATION DIAGNOSTICS FOR QUANTILE MAPPING
# ==============================================================================

#' Compute Diagnostics for Precipitation Quantile Mapping
#'
#' @description
#' Computes diagnostics comparing reference and quantile-mapped precipitation series.
#' Reports moment changes, quantiles, tail metrics, dry/wet-day frequencies, spell lengths,
#' and optional temporal/monthly patterns. When `mean_factor`/`var_factor` are provided
#' together with `month` and `year`, moment and monthly diagnostics are compared against
#' the intended perturbations (scenario-aware null hypothesis).
#'
#' @param precip_ref Numeric vector. Reference (original) precipitation values.
#' @param precip_adj Numeric vector. Adjusted precipitation values (same length as `precip_ref`).
#' @param month Integer vector (optional). Month index (1--12) for each observation.
#' @param year Integer vector (optional). Year index (1..n_years) for each observation.
#' @param mean_factor Numeric matrix (optional). Target mean scaling factors (n_years x 12).
#' @param var_factor Numeric matrix (optional). Target variance scaling factors (n_years x 12).
#' @param wet_thresh Numeric. Threshold for defining wet days (default = 0.1 mm).
#' @param probs Numeric vector. Quantile probabilities to evaluate.
#'
#' @return A list of class `precip_qm_diagnostics` with elements:
#' \itemize{
#'   \item `moments`: moment diagnostics table
#'   \item `quantiles`: quantile comparison table
#'   \item `extremes`: tail metrics table
#'   \item `temporal`: temporal metrics table (optional)
#'   \item `monthly`: month-by-month diagnostics table (optional)
#'   \item `spells`: wet/dry spell diagnostics table
#'   \item `drydays`: wet/dry frequency diagnostics table
#'   \item `summary`: compact summary list
#' }
#'
#' @export
diagnose_precip_qm <- function(
    precip_ref,
    precip_adj,
    month = NULL,
    year = NULL,
    mean_factor = NULL,
    var_factor = NULL,
    wet_thresh = 0.1,
    probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
) {

  # --------------------------------------------------------------------------
  # INPUT VALIDATION
  # --------------------------------------------------------------------------

  if (length(precip_ref) != length(precip_adj)) {
    stop("'precip_ref' and 'precip_adj' must have same length", call. = FALSE)
  }

  n <- length(precip_ref)

  if (!is.null(month) && length(month) != n) {
    stop("'month' must have same length as precipitation series", call. = FALSE)
  }
  if (!is.null(year) && length(year) != n) {
    stop("'year' must have same length as precipitation series", call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # MASKS
  # --------------------------------------------------------------------------

  is_nz_ref  <- precip_ref > 0
  is_nz_adj  <- precip_adj > 0
  is_nz_both <- is_nz_ref & is_nz_adj

  # Note: spell/dryday diagnostics use wet_thresh; moments/quantiles/extremes use nonzero-both mask.

  # --------------------------------------------------------------------------
  # 1) MOMENTS (scenario-aware when factors + month/year are available)
  # --------------------------------------------------------------------------

  moments <- compute_moment_diagnostics(
    precip_ref = precip_ref,
    precip_adj = precip_adj,
    mask = is_nz_both,
    month = month,
    year = year,
    mean_factor = mean_factor,
    var_factor = var_factor
  )

  # --------------------------------------------------------------------------
  # 2) QUANTILES (no intended target unless you define one; reports change)
  # --------------------------------------------------------------------------

  quant_metrics <- compute_quantile_diagnostics(
    precip_ref = precip_ref[is_nz_both],
    precip_adj = precip_adj[is_nz_both],
    probs = probs
  )

  # --------------------------------------------------------------------------
  # 3) EXTREMES (tail ratios)
  # --------------------------------------------------------------------------

  extreme_metrics <- compute_extreme_diagnostics(
    precip_ref = precip_ref[is_nz_both],
    precip_adj = precip_adj[is_nz_both]
  )

  # --------------------------------------------------------------------------
  # 4) TEMPORAL (optional)
  # --------------------------------------------------------------------------

  temporal_metrics <- NULL
  if (!is.null(month) && !is.null(year)) {
    temporal_metrics <- compute_temporal_diagnostics(
      precip_ref = precip_ref,
      precip_adj = precip_adj,
      month = month,
      year = year,
      mask = is_nz_both
    )
  }

  # --------------------------------------------------------------------------
  # 5) MONTHLY (optional; scenario-aware when factors + year are available)
  # --------------------------------------------------------------------------

  monthly_metrics <- NULL
  if (!is.null(month)) {
    monthly_metrics <- compute_monthly_diagnostics(
      precip_ref = precip_ref,
      precip_adj = precip_adj,
      month = month,
      mean_factor = mean_factor,
      var_factor = var_factor,
      year = year
    )
  }

  # --------------------------------------------------------------------------
  # 6) SPELLS (intended preserved unless occurrence is changed)
  # --------------------------------------------------------------------------

  spell_metrics <- compute_spell_diagnostics(
    precip_ref = precip_ref,
    precip_adj = precip_adj,
    wet_thresh = wet_thresh
  )

  # --------------------------------------------------------------------------
  # 7) DRY/WET DAYS (intended preserved for your QM design)
  # --------------------------------------------------------------------------

  dryday_metrics <- compute_dryday_diagnostics(
    precip_ref = precip_ref,
    precip_adj = precip_adj,
    wet_thresh = wet_thresh
  )

  # --------------------------------------------------------------------------
  # 8) SUMMARY (no scoring; scenario-aware where possible)
  # --------------------------------------------------------------------------

  summary_metrics <- compute_summary_metrics(
    moments = moments,
    quant_metrics = quant_metrics,
    extreme_metrics = extreme_metrics,
    spell_metrics = spell_metrics,
    dryday_metrics = dryday_metrics
  )

  result <- list(
    moments = moments,
    quantiles = quant_metrics,
    extremes = extreme_metrics,
    temporal = temporal_metrics,
    monthly = monthly_metrics,
    spells = spell_metrics,
    drydays = dryday_metrics,
    summary = summary_metrics
  )

  class(result) <- c("precip_qm_diagnostics", "list")
  result
}


#' Validate precipitation quantile mapping adjustments
#'
#' @description
#' Convenience wrapper around `diagnose_precip_qm()` with argument names aligned to QM workflows.
#'
#' @param precip_org Numeric vector. Original precipitation values.
#' @param precip_adjusted Numeric vector. Adjusted precipitation values (same length as `precip_org`).
#' @param month Integer vector (optional). Month index (1--12) for each observation.
#' @param year Integer vector (optional). Year index (1..n_years) for each observation.
#' @param mean_factor Numeric matrix (optional). Target mean scaling factors (n_years x 12).
#' @param var_factor Numeric matrix (optional). Target variance scaling factors (n_years x 12).
#' @param wet_thresh Numeric. Threshold for defining wet days (default = 0.1 mm).
#' @param probs Numeric vector. Quantile probabilities to evaluate.
#'
#' @return A `precip_qm_diagnostics` object (see `diagnose_precip_qm()`).
#'
#' @export
validate_quantile_mapping <- function(
    precip_org,
    precip_adjusted,
    month = NULL,
    year = NULL,
    mean_factor = NULL,
    var_factor = NULL,
    wet_thresh = 0.1,
    probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
) {

  if (length(precip_org) != length(precip_adjusted)) {
    stop("'precip_org' and 'precip_adjusted' must have same length", call. = FALSE)
  }

  diagnose_precip_qm(
    precip_ref = precip_org,
    precip_adj = precip_adjusted,
    month = month,
    year = year,
    mean_factor = mean_factor,
    var_factor = var_factor,
    wet_thresh = wet_thresh,
    probs = probs
  )
}

# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Compute intended scaling ratios over the evaluated subset (time-weighted)
#' @keywords internal
.compute_intended_ratios <- function(mask, month, year, mean_factor, var_factor) {

  out <- list(mean_ratio = NA_real_, var_ratio = NA_real_)

  if (is.null(month) || is.null(year)) return(out)
  if (is.null(mean_factor) || is.null(var_factor)) return(out)
  if (!is.matrix(mean_factor) || !is.matrix(var_factor)) return(out)

  idx <- which(mask)
  if (length(idx) == 0) return(out)

  # Expected factor per observation (year-month)
  f_mean <- mean_factor[cbind(year[idx], month[idx])]
  f_var  <- var_factor[cbind(year[idx], month[idx])]

  out$mean_ratio <- mean(f_mean, na.rm = TRUE)
  out$var_ratio  <- mean(f_var,  na.rm = TRUE)

  out
}

# ==============================================================================
# DIAGNOSTIC TABLE BUILDERS
# ==============================================================================

#' Compute moment diagnostics (scenario-aware when factors + month/year are available)
#' @keywords internal
compute_moment_diagnostics <- function(precip_ref, precip_adj, mask,
                                       month = NULL, year = NULL,
                                       mean_factor = NULL, var_factor = NULL,
                                       tol_mean = 0.05, tol_var = 0.10) {

  ref <- precip_ref[mask]
  adj <- precip_adj[mask]

  overall <- data.frame(
    metric = c("mean", "sd", "variance", "cv", "skewness", "kurtosis"),
    original = c(
      mean(ref),
      sd(ref),
      stats::var(ref),
      sd(ref) / mean(ref),
      compute_skewness(ref),
      compute_kurtosis(ref)
    ),
    adjusted = c(
      mean(adj),
      sd(adj),
      stats::var(adj),
      sd(adj) / mean(adj),
      compute_skewness(adj),
      compute_kurtosis(adj)
    ),
    stringsAsFactors = FALSE
  )

  overall$ratio <- overall$adjusted / overall$original
  overall$diff  <- overall$adjusted - overall$original
  overall$pct_change <- (overall$ratio - 1) * 100

  intended <- .compute_intended_ratios(
    mask = mask,
    month = month,
    year = year,
    mean_factor = mean_factor,
    var_factor = var_factor
  )

  overall$intended_ratio <- NA_real_
  overall$ratio_error <- NA_real_
  overall$pct_error <- NA_real_
  overall$within_tolerance <- NA

  # mean target
  i_mean <- which(overall$metric == "mean")
  if (is.finite(intended$mean_ratio)) {
    overall$intended_ratio[i_mean] <- intended$mean_ratio
    overall$ratio_error[i_mean] <- overall$ratio[i_mean] - intended$mean_ratio
    overall$pct_error[i_mean] <- 100 * overall$ratio_error[i_mean] / intended$mean_ratio
    overall$within_tolerance[i_mean] <- abs(overall$ratio_error[i_mean]) <= tol_mean
  }

  # variance target
  i_var <- which(overall$metric == "variance")
  if (is.finite(intended$var_ratio)) {
    overall$intended_ratio[i_var] <- intended$var_ratio
    overall$ratio_error[i_var] <- overall$ratio[i_var] - intended$var_ratio
    overall$pct_error[i_var] <- 100 * overall$ratio_error[i_var] / intended$var_ratio
    overall$within_tolerance[i_var] <- abs(overall$ratio_error[i_var]) <= tol_var
  }

  # sd target derived from variance target
  i_sd <- which(overall$metric == "sd")
  if (is.finite(intended$var_ratio)) {
    sd_intended <- sqrt(intended$var_ratio)
    overall$intended_ratio[i_sd] <- sd_intended
    overall$ratio_error[i_sd] <- overall$ratio[i_sd] - sd_intended
    overall$pct_error[i_sd] <- 100 * overall$ratio_error[i_sd] / sd_intended
    overall$within_tolerance[i_sd] <- abs(overall$ratio_error[i_sd]) <= sqrt(tol_var)
  }

  overall
}

#' Compute quantile diagnostics (reports change; no scenario-aware target by default)
#' @keywords internal
compute_quantile_diagnostics <- function(precip_ref, precip_adj, probs) {

  q_ref <- stats::quantile(precip_ref, probs = probs, na.rm = TRUE)
  q_adj <- stats::quantile(precip_adj, probs = probs, na.rm = TRUE)

  result <- data.frame(
    prob = probs,
    original = as.numeric(q_ref),
    adjusted = as.numeric(q_adj),
    stringsAsFactors = FALSE
  )

  result$ratio <- result$adjusted / result$original
  result$diff  <- result$adjusted - result$original
  result$pct_change <- (result$ratio - 1) * 100
  result$abs_error <- abs(result$diff)

  result
}

#' Compute extreme-value diagnostics (tail ratios)
#' @keywords internal
compute_extreme_diagnostics <- function(precip_ref, precip_adj) {

  thresholds <- c(0.90, 0.95, 0.99, 0.999)

  metrics <- data.frame(
    threshold = thresholds,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(thresholds)) {
    p <- thresholds[i]

    q_ref <- stats::quantile(precip_ref, p, na.rm = TRUE)
    tail_ref <- precip_ref[precip_ref >= q_ref]

    q_adj <- stats::quantile(precip_adj, p, na.rm = TRUE)
    tail_adj <- precip_adj[precip_adj >= q_adj]

    metrics$original_threshold[i] <- q_ref
    metrics$adjusted_threshold[i] <- q_adj
    metrics$original_mean[i] <- mean(tail_ref, na.rm = TRUE)
    metrics$adjusted_mean[i] <- mean(tail_adj, na.rm = TRUE)
    metrics$original_max[i] <- max(tail_ref, na.rm = TRUE)
    metrics$adjusted_max[i] <- max(tail_adj, na.rm = TRUE)
  }

  metrics$threshold_ratio <- metrics$adjusted_threshold / metrics$original_threshold
  metrics$mean_ratio <- metrics$adjusted_mean / metrics$original_mean
  metrics$max_ratio <- metrics$adjusted_max / metrics$original_max

  metrics
}

#' Compute temporal diagnostics (no subjective categories)
#' @keywords internal
compute_temporal_diagnostics <- function(precip_ref, precip_adj, month, year, mask) {

  annual_ref <- tapply(precip_ref[mask], year[mask], mean, na.rm = TRUE)
  annual_adj <- tapply(precip_adj[mask], year[mask], mean, na.rm = TRUE)

  monthly_ref <- tapply(precip_ref[mask], month[mask], mean, na.rm = TRUE)
  monthly_adj <- tapply(precip_adj[mask], month[mask], mean, na.rm = TRUE)

  cor_annual  <- stats::cor(annual_ref,  annual_adj,  use = "complete.obs")
  cor_monthly <- stats::cor(monthly_ref, monthly_adj, use = "complete.obs")

  acf_ref <- stats::acf(precip_ref[mask], lag.max = 1, plot = FALSE)$acf[2]
  acf_adj <- stats::acf(precip_adj[mask], lag.max = 1, plot = FALSE)$acf[2]

  data.frame(
    metric = c(
      "annual_correlation",
      "monthly_correlation",
      "lag1_autocorr_ref",
      "lag1_autocorr_adj",
      "autocorr_abs_diff"
    ),
    value = c(
      cor_annual,
      cor_monthly,
      acf_ref,
      acf_adj,
      abs(acf_adj - acf_ref)
    ),
    stringsAsFactors = FALSE
  )
}

#' Compute month-by-month diagnostics (scenario-aware when factors + year are available)
#' @keywords internal
compute_monthly_diagnostics <- function(precip_ref, precip_adj, month,
                                        mean_factor = NULL, var_factor = NULL,
                                        year = NULL,
                                        tol_mean = 0.05, tol_var = 0.10) {

  monthly <- data.frame(
    month = 1:12,
    stringsAsFactors = FALSE
  )

  for (m in 1:12) {
    idx <- month == m & precip_ref > 0 & precip_adj > 0

    if (sum(idx) > 0) {
      monthly$original_mean[m] <- mean(precip_ref[idx], na.rm = TRUE)
      monthly$adjusted_mean[m] <- mean(precip_adj[idx], na.rm = TRUE)
      monthly$original_var[m]  <- stats::var(precip_ref[idx], na.rm = TRUE)
      monthly$adjusted_var[m]  <- stats::var(precip_adj[idx], na.rm = TRUE)
      monthly$n_obs[m] <- sum(idx)
    } else {
      monthly$original_mean[m] <- NA_real_
      monthly$adjusted_mean[m] <- NA_real_
      monthly$original_var[m]  <- NA_real_
      monthly$adjusted_var[m]  <- NA_real_
      monthly$n_obs[m] <- 0
    }
  }

  monthly$mean_ratio <- monthly$adjusted_mean / monthly$original_mean
  monthly$var_ratio  <- monthly$adjusted_var  / monthly$original_var

  if (!is.null(mean_factor) && !is.null(var_factor) && !is.null(year)) {

    monthly$target_mean_ratio <- colMeans(mean_factor, na.rm = TRUE)
    monthly$target_var_ratio  <- colMeans(var_factor,  na.rm = TRUE)

    monthly$mean_error <- monthly$mean_ratio - monthly$target_mean_ratio
    monthly$var_error  <- monthly$var_ratio  - monthly$target_var_ratio

    monthly$abs_mean_error <- abs(monthly$mean_error)
    monthly$abs_var_error  <- abs(monthly$var_error)

    monthly$within_tol_mean <- monthly$abs_mean_error <= tol_mean
    monthly$within_tol_var  <- monthly$abs_var_error  <= tol_var
  }

  monthly
}

#' Compute wet and dry spell-length diagnostics (intended preserved)
#' @keywords internal
#' Compute wet and dry spell-length diagnostics (intended preserved)
#' @keywords internal
compute_spell_diagnostics <- function(precip_ref, precip_adj, wet_thresh,
                                      tol_mean_ratio = 0.20) {

  .safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x, na.rm = TRUE)
  }

  .safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
  }

  .safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    max(x, na.rm = TRUE)
  }

  dry_ref <- compute_spell_lengths(precip_ref, wet_thresh, below = TRUE)
  dry_adj <- compute_spell_lengths(precip_adj, wet_thresh, below = TRUE)

  wet_ref <- compute_spell_lengths(precip_ref, wet_thresh, below = FALSE)
  wet_adj <- compute_spell_lengths(precip_adj, wet_thresh, below = FALSE)

  result <- data.frame(
    spell_type = c("dry", "wet"),
    original_mean = c(.safe_mean(dry_ref), .safe_mean(wet_ref)),
    adjusted_mean = c(.safe_mean(dry_adj), .safe_mean(wet_adj)),
    original_max  = c(.safe_max(dry_ref),  .safe_max(wet_ref)),
    adjusted_max  = c(.safe_max(dry_adj),  .safe_max(wet_adj)),
    original_sd   = c(.safe_sd(dry_ref),   .safe_sd(wet_ref)),
    adjusted_sd   = c(.safe_sd(dry_adj),   .safe_sd(wet_adj)),
    stringsAsFactors = FALSE
  )

  result$mean_ratio <- result$adjusted_mean / result$original_mean
  result$max_ratio  <- result$adjusted_max  / result$original_max
  result$sd_ratio   <- result$adjusted_sd   / result$original_sd

  result$intended_mean_ratio <- 1.0
  result$mean_ratio_error <- result$mean_ratio - result$intended_mean_ratio

  # If mean_ratio is NA (e.g., no spells), within_tolerance should be NA, not FALSE
  result$within_tolerance <- ifelse(
    is.finite(result$mean_ratio_error),
    abs(result$mean_ratio_error) <= tol_mean_ratio,
    NA
  )

  result
}


#' Compute dry-day and wet-day frequency diagnostics (intended preserved)
#' @keywords internal
compute_dryday_diagnostics <- function(precip_ref, precip_adj, wet_thresh,
                                       tol_freq_abs = 0.01) {

  n_total <- length(precip_ref)

  dry_ref <- sum(precip_ref < wet_thresh)
  dry_adj <- sum(precip_adj < wet_thresh)

  wet_ref <- sum(precip_ref >= wet_thresh)
  wet_adj <- sum(precip_adj >= wet_thresh)

  result <- data.frame(
    category = c("dry_days", "wet_days", "dry_freq", "wet_freq"),
    original = c(dry_ref, wet_ref, dry_ref / n_total, wet_ref / n_total),
    adjusted = c(dry_adj, wet_adj, dry_adj / n_total, wet_adj / n_total),
    stringsAsFactors = FALSE
  )

  result$diff <- result$adjusted - result$original
  result$intended_diff <- 0

  # Ratios can be undefined when original == 0; keep but do not enforce
  result$ratio <- suppressWarnings(result$adjusted / result$original)
  result$pct_change <- suppressWarnings((result$ratio - 1) * 100)

  result$within_tolerance <- c(
    result$diff[result$category == "dry_days"] == 0,
    result$diff[result$category == "wet_days"] == 0,
    abs(result$diff[result$category == "dry_freq"]) <= tol_freq_abs,
    abs(result$diff[result$category == "wet_freq"]) <= tol_freq_abs
  )

  result
}

#' Compute summary metrics (no qualitative grades; no scoring)
#' @keywords internal
compute_summary_metrics <- function(moments, quant_metrics, extreme_metrics,
                                    spell_metrics, dryday_metrics) {

  out <- list()

  mean_row <- moments[moments$metric == "mean", , drop = FALSE]
  var_row  <- moments[moments$metric == "variance", , drop = FALSE]

  out$moments <- list(
    mean_ratio_observed   = mean_row$ratio,
    mean_ratio_intended   = mean_row$intended_ratio,
    mean_ratio_error      = mean_row$ratio_error,
    mean_within_tolerance = mean_row$within_tolerance,

    var_ratio_observed    = var_row$ratio,
    var_ratio_intended    = var_row$intended_ratio,
    var_ratio_error       = var_row$ratio_error,
    var_within_tolerance  = var_row$within_tolerance
  )

  out$quantiles <- list(
    max_abs_pct_change = max(abs(quant_metrics$pct_change), na.rm = TRUE),
    max_abs_ratio_minus_1 = max(abs(quant_metrics$ratio - 1), na.rm = TRUE)
  )

  out$extremes <- list(
    mean_threshold_ratio = mean(extreme_metrics$threshold_ratio, na.rm = TRUE),
    mean_tail_mean_ratio = mean(extreme_metrics$mean_ratio, na.rm = TRUE),
    mean_max_ratio       = mean(extreme_metrics$max_ratio, na.rm = TRUE)
  )

  dd_dry <- dryday_metrics[dryday_metrics$category == "dry_days", , drop = FALSE]
  dd_wet <- dryday_metrics[dryday_metrics$category == "wet_days", , drop = FALSE]

  out$drydays <- list(
    dry_days_diff = dd_dry$diff,
    wet_days_diff = dd_wet$diff,
    dry_days_within_tolerance = dd_dry$within_tolerance,
    wet_days_within_tolerance = dd_wet$within_tolerance
  )

  out$spells <- list(
    dry_spell_mean_ratio = spell_metrics$mean_ratio[spell_metrics$spell_type == "dry"],
    wet_spell_mean_ratio = spell_metrics$mean_ratio[spell_metrics$spell_type == "wet"],
    dry_spell_within_tolerance = spell_metrics$within_tolerance[spell_metrics$spell_type == "dry"],
    wet_spell_within_tolerance = spell_metrics$within_tolerance[spell_metrics$spell_type == "wet"]
  )

  out
}

# ==============================================================================
# PRINT HELPER FUNCTIONS
# ==============================================================================

#' Build moments table for printing
#' @keywords internal
.make_moments_table <- function(moments) {
  keep <- moments$metric %in% c("mean", "variance")
  m <- moments[keep, , drop = FALSE]

  target_pct <- 100 * (m$intended_ratio - 1)
  achieved_pct <- 100 * (m$ratio - 1)
  error_pp <- achieved_pct - target_pct

  data.frame(
    Statistic = c(mean = "Mean", variance = "Variance")[m$metric],
    Original = .format_num(m$original, 2),
    Adjusted = .format_num(m$adjusted, 2),
    Target = .format_pct(target_pct, 0),
    Achieved = .format_pct(achieved_pct, 0),
    Error = .format_pct(error_pp, 0),
    stringsAsFactors = FALSE
  )
}

#' Build dry/wet day count table for printing
#' @keywords internal
.make_drywet_table <- function(drydays) {
  dd <- drydays[drydays$category %in% c("dry_days", "wet_days"), , drop = FALSE]
  lab <- c(dry_days = "Dry days", wet_days = "Wet days")

  data.frame(
    Category = lab[dd$category],
    Original = as.integer(dd$original),
    Adjusted = as.integer(dd$adjusted),
    Difference = as.integer(dd$diff),
    stringsAsFactors = FALSE
  )
}

#' Build quantiles table for printing
#' @keywords internal
.make_quantiles_table <- function(quantiles) {
  keep_p <- c(0.50, 0.90, 0.95, 0.99)
  q <- quantiles[quantiles$prob %in% keep_p, , drop = FALSE]

  lab <- c(
    "0.5" = "P50 (Median)",
    "0.9" = "P90",
    "0.95" = "P95",
    "0.99" = "P99"
  )

  pct <- 100 * (q$ratio - 1)

  data.frame(
    Percentile = lab[as.character(q$prob)],
    Original = .format_num(q$original, 2),
    Adjusted = .format_num(q$adjusted, 2),
    Change = .format_pct(pct, 0),
    stringsAsFactors = FALSE
  )
}

#' Build extremes table for printing
#' @keywords internal
.make_extremes_table <- function(extremes) {
  if (is.null(extremes) || !is.data.frame(extremes) || nrow(extremes) == 0) {
    return(NULL)
  }
  if (!("threshold" %in% names(extremes))) {
    return(NULL)
  }

  ex <- extremes

  # Build readable labels from threshold values
  p <- ex$threshold
  pct <- round(100 * p, 1)
  top <- round(100 * (1 - p), 3)

  fmt_pctile <- function(v) {
    ifelse(abs(v - round(v)) < 1e-8, sprintf("%.0f", round(v)), sprintf("%.1f", v))
  }
  fmt_top <- function(v) {
    ifelse(v < 1, sprintf("%.1f%%", v), sprintf("%.0f%%", v))
  }

  exceedance_label <- paste0("Top ", fmt_top(top), " (P", fmt_pctile(pct), ")")

  data.frame(
    Tail = exceedance_label,
    `Threshold change` = .format_pct(100 * (ex$threshold_ratio - 1), 0),
    `Mean change` = .format_pct(100 * (ex$mean_ratio - 1), 0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Build spell lengths table for printing
#' @keywords internal
.make_spells_table <- function(spells) {
  lab <- c(dry = "Dry spells", wet = "Wet spells")
  pct <- 100 * (spells$mean_ratio - 1)

  data.frame(
    Type = lab[spells$spell_type],
    `Original (days)` = .format_num(spells$original_mean, 1),
    `Adjusted (days)` = .format_num(spells$adjusted_mean, 1),
    Change = .format_pct(pct, 0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Build monthly table for printing
#' @keywords internal
.make_monthly_table <- function(monthly) {
  if (is.null(monthly)) return(NULL)
  if (!all(c("target_mean_ratio", "mean_ratio") %in% names(monthly))) return(NULL)

  target <- 100 * (monthly$target_mean_ratio - 1)
  achieved <- 100 * (monthly$mean_ratio - 1)
  error <- achieved - target

  data.frame(
    Month = month.abb[monthly$month],
    Target = .format_pct(target, 0),
    Achieved = .format_pct(achieved, 0),
    Error = .format_pct(error, 0),
    stringsAsFactors = FALSE
  )
}


#' Print method for precipitation QM diagnostics
#'
#' @description
#' Prints a formatted summary of precipitation quantile mapping diagnostics
#' as multiple tables showing key metrics.
#'
#' @param x An object of class `precip_qm_diagnostics`.
#' @param ... Not used. Included for S3 method consistency.
#'
#' @return The input `x`, invisibly.
#'
#' @export
print.precip_qm_diagnostics <- function(x, ...) {


  cat("\n=== Precipitation QM Diagnostics ===\n\n")

  # Table 1: Overall moment changes
  cat("Mean and variance changes:\n")
  print(.make_moments_table(x$moments), row.names = FALSE)
  cat("\n")

  # Table 2: Dry/wet day counts
  cat("Wet/dry day counts:\n")
  print(.make_drywet_table(x$drydays), row.names = FALSE)
  cat("\n")

  # Table 3: Quantiles
  cat("Percentiles (wet days):\n")
  print(.make_quantiles_table(x$quantiles), row.names = FALSE)
  cat("\n")

  # Table 4: Extremes
  ext_tab <- .make_extremes_table(x$extremes)
  if (!is.null(ext_tab)) {
    cat("Extreme tails:\n")
    print(ext_tab, row.names = FALSE)
    cat("\n")
  }

  # Table 5: Spell lengths
  cat("Spell lengths:\n")
  print(.make_spells_table(x$spells), row.names = FALSE)
  cat("\n")

  # Table 6: Monthly (if available)
  mtab <- .make_monthly_table(x$monthly)
  if (!is.null(mtab)) {
    cat("Monthly mean changes:\n")
    print(mtab, row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}


#' Summarize precipitation QM diagnostics
#'
#' @description
#' Prints a compact numeric summary of the key QM diagnostic metrics.
#'
#' @param object A `precip_qm_diagnostics` object.
#' @param ... Unused.
#'
#' @return The input `object`, invisibly.
#'
#' @export
summary.precip_qm_diagnostics <- function(object, ...) {

  cat("\n=== Precipitation QM Summary ===\n\n")

  s <- object$summary

  # Moments (scenario-aware)
  cat("Targets vs achieved:\n")
  if (is.finite(s$moments$mean_ratio_observed)) {
    cat(sprintf("  Mean:     achieved %.3f", s$moments$mean_ratio_observed))
    if (is.finite(s$moments$mean_ratio_intended)) {
      cat(sprintf(" | target %.3f | error %+.3f\n",
                  s$moments$mean_ratio_intended,
                  s$moments$mean_ratio_error))
    } else {
      cat("\n")
    }
  }

  if (is.finite(s$moments$var_ratio_observed)) {
    cat(sprintf("  Variance: achieved %.3f", s$moments$var_ratio_observed))
    if (is.finite(s$moments$var_ratio_intended)) {
      cat(sprintf(" | target %.3f | error %+.3f\n",
                  s$moments$var_ratio_intended,
                  s$moments$var_ratio_error))
    } else {
      cat("\n")
    }
  }
  cat("\n")

  # Dry/wet days
  cat("Wet/dry day preservation:\n")
  cat(sprintf("  Dry days change: %d\n", as.integer(s$drydays$dry_days_diff)))
  cat(sprintf("  Wet days change: %d\n", as.integer(s$drydays$wet_days_diff)))
  cat("\n")

  # Quantiles
  cat("Quantile changes:\n")
  cat(sprintf("  Max absolute change: %.1f%%\n", s$quantiles$max_abs_pct_change))
  cat("\n")

  # Extremes
  cat("Extreme tail ratios (mean across thresholds):\n")
  cat(sprintf("  Threshold: %.3f\n", s$extremes$mean_threshold_ratio))
  cat(sprintf("  Tail mean: %.3f\n", s$extremes$mean_tail_mean_ratio))
  cat("\n")

  invisible(object)
}
