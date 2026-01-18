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
#' @param prcp_ref Numeric vector. Reference (original) precipitation values.
#' @param prcp_adj Numeric vector. Adjusted precipitation values (same length as `prcp_ref`).
#' @param month Integer vector (optional). Month index (1--12) for each observation.
#' @param year Integer vector (optional). Year index (1..n_years) for each observation.
#' @param mean_factor Numeric matrix (optional). Target mean scaling factors (n_years x 12).
#' @param var_factor Numeric matrix (optional). Target variance scaling factors (n_years x 12).
#' @param wet_thresh Numeric. Threshold for defining wet days (default = 0.1 mm).
#' @param probs Numeric vector. Quantile probabilities to evaluate.
#'
#' @return A list of class `prcp_qm_diagnostics` with elements:
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
diagnose_prcp_qm <- function(
    prcp_ref,
    prcp_adj,
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

  if (length(prcp_ref) != length(prcp_adj)) {
    stop("'prcp_ref' and 'prcp_adj' must have same length", call. = FALSE)
  }

  n <- length(prcp_ref)

  if (!is.null(month) && length(month) != n) {
    stop("'month' must have same length as precipitation series", call. = FALSE)
  }
  if (!is.null(year) && length(year) != n) {
    stop("'year' must have same length as precipitation series", call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # MASKS
  # --------------------------------------------------------------------------

  is_nz_ref  <- prcp_ref > 0
  is_nz_adj  <- prcp_adj > 0
  is_nz_both <- is_nz_ref & is_nz_adj

  # Note: spell/dryday diagnostics use wet_thresh; moments/quantiles/extremes use nonzero-both mask.

  # --------------------------------------------------------------------------
  # 1) MOMENTS (scenario-aware when factors + month/year are available)
  # --------------------------------------------------------------------------

  moments <- compute_moment_diagnostics(
    prcp_ref = prcp_ref,
    prcp_adj = prcp_adj,
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
    prcp_ref = prcp_ref[is_nz_both],
    prcp_adj = prcp_adj[is_nz_both],
    probs = probs
  )

  # --------------------------------------------------------------------------
  # 3) EXTREMES (tail ratios)
  # --------------------------------------------------------------------------

  extreme_metrics <- compute_extreme_diagnostics(
    prcp_ref = prcp_ref[is_nz_both],
    prcp_adj = prcp_adj[is_nz_both]
  )

  # --------------------------------------------------------------------------
  # 4) TEMPORAL (optional)
  # --------------------------------------------------------------------------

  temporal_metrics <- NULL
  if (!is.null(month) && !is.null(year)) {
    temporal_metrics <- compute_temporal_diagnostics(
      prcp_ref = prcp_ref,
      prcp_adj = prcp_adj,
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
      prcp_ref = prcp_ref,
      prcp_adj = prcp_adj,
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
    prcp_ref = prcp_ref,
    prcp_adj = prcp_adj,
    wet_thresh = wet_thresh
  )

  # --------------------------------------------------------------------------
  # 7) DRY/WET DAYS (intended preserved for your QM design)
  # --------------------------------------------------------------------------

  dryday_metrics <- compute_dryday_diagnostics(
    prcp_ref = prcp_ref,
    prcp_adj = prcp_adj,
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

  class(result) <- c("prcp_qm_diagnostics", "list")
  result
}


#' Validate precipitation quantile mapping adjustments
#'
#' @description
#' Convenience wrapper around `diagnose_prcp_qm()` with argument names aligned to QM workflows.
#'
#' @param prcp_org Numeric vector. Original precipitation values.
#' @param prcp_adjusted Numeric vector. Adjusted precipitation values (same length as `prcp_org`).
#' @param month Integer vector (optional). Month index (1--12) for each observation.
#' @param year Integer vector (optional). Year index (1..n_years) for each observation.
#' @param mean_factor Numeric matrix (optional). Target mean scaling factors (n_years x 12).
#' @param var_factor Numeric matrix (optional). Target variance scaling factors (n_years x 12).
#' @param wet_thresh Numeric. Threshold for defining wet days (default = 0.1 mm).
#' @param probs Numeric vector. Quantile probabilities to evaluate.
#'
#' @return A `prcp_qm_diagnostics` object (see `diagnose_prcp_qm()`).
#'
#' @export
validate_quantile_mapping <- function(
    prcp_org,
    prcp_adjusted,
    month = NULL,
    year = NULL,
    mean_factor = NULL,
    var_factor = NULL,
    wet_thresh = 0.1,
    probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
) {

  if (length(prcp_org) != length(prcp_adjusted)) {
    stop("'prcp_org' and 'prcp_adjusted' must have same length", call. = FALSE)
  }

  diagnose_prcp_qm(
    prcp_ref = prcp_org,
    prcp_adj = prcp_adjusted,
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
compute_moment_diagnostics <- function(prcp_ref, prcp_adj, mask,
                                       month = NULL, year = NULL,
                                       mean_factor = NULL, var_factor = NULL,
                                       tol_mean = 0.05, tol_var = 0.10) {

  ref <- prcp_ref[mask]
  adj <- prcp_adj[mask]

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
compute_quantile_diagnostics <- function(prcp_ref, prcp_adj, probs) {

  q_ref <- stats::quantile(prcp_ref, probs = probs, na.rm = TRUE)
  q_adj <- stats::quantile(prcp_adj, probs = probs, na.rm = TRUE)

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
compute_extreme_diagnostics <- function(prcp_ref, prcp_adj) {

  thresholds <- c(0.90, 0.95, 0.99, 0.999)

  metrics <- data.frame(
    threshold = thresholds,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(thresholds)) {
    p <- thresholds[i]

    q_ref <- stats::quantile(prcp_ref, p, na.rm = TRUE)
    tail_ref <- prcp_ref[prcp_ref >= q_ref]

    q_adj <- stats::quantile(prcp_adj, p, na.rm = TRUE)
    tail_adj <- prcp_adj[prcp_adj >= q_adj]

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
compute_temporal_diagnostics <- function(prcp_ref, prcp_adj, month, year, mask) {

  annual_ref <- tapply(prcp_ref[mask], year[mask], mean, na.rm = TRUE)
  annual_adj <- tapply(prcp_adj[mask], year[mask], mean, na.rm = TRUE)

  monthly_ref <- tapply(prcp_ref[mask], month[mask], mean, na.rm = TRUE)
  monthly_adj <- tapply(prcp_adj[mask], month[mask], mean, na.rm = TRUE)

  cor_annual  <- stats::cor(annual_ref,  annual_adj,  use = "complete.obs")
  cor_monthly <- stats::cor(monthly_ref, monthly_adj, use = "complete.obs")

  acf_ref <- stats::acf(prcp_ref[mask], lag.max = 1, plot = FALSE)$acf[2]
  acf_adj <- stats::acf(prcp_adj[mask], lag.max = 1, plot = FALSE)$acf[2]

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
compute_monthly_diagnostics <- function(prcp_ref, prcp_adj, month,
                                        mean_factor = NULL, var_factor = NULL,
                                        year = NULL,
                                        tol_mean = 0.05, tol_var = 0.10) {

  monthly <- data.frame(
    month = 1:12,
    stringsAsFactors = FALSE
  )

  for (m in 1:12) {
    idx <- month == m & prcp_ref > 0 & prcp_adj > 0

    if (sum(idx) > 0) {
      monthly$original_mean[m] <- mean(prcp_ref[idx], na.rm = TRUE)
      monthly$adjusted_mean[m] <- mean(prcp_adj[idx], na.rm = TRUE)
      monthly$original_var[m]  <- stats::var(prcp_ref[idx], na.rm = TRUE)
      monthly$adjusted_var[m]  <- stats::var(prcp_adj[idx], na.rm = TRUE)
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
compute_spell_diagnostics <- function(prcp_ref, prcp_adj, wet_thresh,
                                      tol_mean_ratio = 0.20) {

  dry_ref <- compute_spell_lengths(prcp_ref, wet_thresh, below = TRUE)
  dry_adj <- compute_spell_lengths(prcp_adj, wet_thresh, below = TRUE)

  wet_ref <- compute_spell_lengths(prcp_ref, wet_thresh, below = FALSE)
  wet_adj <- compute_spell_lengths(prcp_adj, wet_thresh, below = FALSE)

  result <- data.frame(
    spell_type = c("dry", "wet"),
    original_mean = c(mean(dry_ref), mean(wet_ref)),
    adjusted_mean = c(mean(dry_adj), mean(wet_adj)),
    original_max = c(max(dry_ref), max(wet_ref)),
    adjusted_max = c(max(dry_adj), max(wet_adj)),
    original_sd = c(sd(dry_ref), sd(wet_ref)),
    adjusted_sd = c(sd(dry_adj), sd(wet_adj)),
    stringsAsFactors = FALSE
  )

  result$mean_ratio <- result$adjusted_mean / result$original_mean
  result$max_ratio  <- result$adjusted_max  / result$original_max
  result$sd_ratio   <- result$adjusted_sd   / result$original_sd

  # Intended preserve unless you intentionally modify occurrence
  result$intended_mean_ratio <- 1.0
  result$mean_ratio_error <- result$mean_ratio - result$intended_mean_ratio
  result$within_tolerance <- abs(result$mean_ratio_error) <= tol_mean_ratio

  result
}

#' Compute dry-day and wet-day frequency diagnostics (intended preserved)
#' @keywords internal
compute_dryday_diagnostics <- function(prcp_ref, prcp_adj, wet_thresh,
                                       tol_freq_abs = 0.01) {

  n_total <- length(prcp_ref)

  dry_ref <- sum(prcp_ref < wet_thresh)
  dry_adj <- sum(prcp_adj < wet_thresh)

  wet_ref <- sum(prcp_ref >= wet_thresh)
  wet_adj <- sum(prcp_adj >= wet_thresh)

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
# PRINT, PLOT, SUMMARY METHODS
# ==============================================================================

.make_moments_table <- function(moments) {
  keep <- moments$metric %in% c("mean", "variance", "sd")
  m <- moments[keep, , drop = FALSE]

  intended_pct <- 100 * (m$intended_ratio - 1)
  observed_pct <- 100 * (m$ratio - 1)
  diff_pp <- observed_pct - intended_pct

  data.frame(
    Metric = c(mean = "Mean", variance = "Variance", sd = "Std. dev.")[m$metric],
    Original = .format_num(m$original, 2),
    Adjusted = .format_num(m$adjusted, 2),
    `Intended change` = .format_pct(intended_pct, 0),
    `Observed change` = .format_pct(observed_pct, 0),
    `Difference (pp)` = .format_pct(diff_pp, 0),
    stringsAsFactors = FALSE
  )
}

.make_drywet_table <- function(drydays) {
  dd <- drydays[drydays$category %in% c("dry_days", "wet_days"), , drop = FALSE]
  lab <- c(dry_days = "Dry days", wet_days = "Wet days")

  data.frame(
    Category = lab[dd$category],
    Original = as.integer(dd$original),
    Adjusted = as.integer(dd$adjusted),
    Change = as.integer(dd$diff),
    stringsAsFactors = FALSE
  )
}

.make_quantiles_table <- function(quantiles) {
  keep_p <- c(0.50, 0.90, 0.95, 0.99)
  q <- quantiles[quantiles$prob %in% keep_p, , drop = FALSE]

  lab <- c(
    "0.5" = "Percentile 50 (Median)",
    "0.9" = "Percentile 90",
    "0.95" = "Percentile 95",
    "0.99" = "Percentile 99"
  )

  pct <- 100 * (q$ratio - 1)

  data.frame(
    Percentile = lab[as.character(q$prob)],
    Original = .format_num(q$original, 2),
    Adjusted = .format_num(q$adjusted, 2),
    `Change` = .format_pct(pct, 0),
    stringsAsFactors = FALSE
  )
}

.make_spells_table <- function(spells) {
  lab <- c(dry = "Dry spell length", wet = "Wet spell length")
  pct <- 100 * (spells$mean_ratio - 1)

  data.frame(
    `Spell type` = lab[spells$spell_type],
    `Original mean (days)` = .format_num(spells$original_mean, 2),
    `Adjusted mean (days)` = .format_num(spells$adjusted_mean, 2),
    Change = .format_pct(pct, 0),
    stringsAsFactors = FALSE
  )
}

.make_monthly_table <- function(monthly) {
  if (is.null(monthly)) return(NULL)
  if (!all(c("target_mean_ratio", "mean_ratio") %in% names(monthly))) return(NULL)

  intended <- 100 * (monthly$target_mean_ratio - 1)
  observed <- 100 * (monthly$mean_ratio - 1)
  diff_pp <- observed - intended

  data.frame(
    Month = month.abb[monthly$month],
    `Intended change` = .format_pct(intended, 0),
    `Observed change` = .format_pct(observed, 0),
    `Difference (pp)` = .format_pct(diff_pp, 0),
    stringsAsFactors = FALSE
  )
}


#' Print method for precipitation QM diagnostics
#'
#' @param x An object of class \code{"prcp_qm_diagnostics"}.
#' @param ... Not used. Included for S3 method consistency.
#'
#' @return The input \code{x}, invisibly.
#' @export
print.prcp_qm_diagnostics <- function(x, ...) {

  cat("\n=== Precipitation QM Diagnostics ===\n\n")

  cat("Overall changes vs targets:\n")
  print(.make_moments_table(x$moments), row.names = FALSE)
  cat("\n")

  cat("Dry and wet day counts:\n")
  print(.make_drywet_table(x$drydays), row.names = FALSE)
  cat("\n")

  cat("Selected percentiles (wet days only):\n")
  print(.make_quantiles_table(x$quantiles), row.names = FALSE)
  cat("\n")

  cat("Extreme precipitation changes (distributional comparison):\n")
  ex <- x$extremes

  ex$exceedance_label <- c(
    "Top 10% (Percentile 90)",
    "Top 5% (Percentile 95)",
    "Top 1% (Percentile 99)",
    "Top 0.1% (Percentile 99.9)"
  )

  ext_print <- data.frame(
    `Exceedance level` = ex$exceedance_label,
    `Extreme cutoff change` = .format_pct(100 * (ex$threshold_ratio - 1), 0),
    `Average extreme size change` = .format_pct(100 * (ex$mean_ratio - 1), 0),
    `Maximum event change` = .format_pct(100 * (ex$max_ratio - 1), 0),
    stringsAsFactors = FALSE
  )
  print(ext_print, row.names = FALSE)
  cat("\n")

  cat("Spell lengths:\n")
  print(.make_spells_table(x$spells), row.names = FALSE)
  cat("\n")

  mtab <- .make_monthly_table(x$monthly)
  if (!is.null(mtab)) {
    cat("Monthly mean change vs target:\n")
    print(mtab, row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}


#' Print precipitation QM diagnostics
#' @export
print.prcp_qm_diagnostics <- function(x, ...) {

  cat("\n=== Precipitation QM Diagnostics ===\n\n")

  cat("Overall changes vs targets:\n")
  print(.make_moments_table(x$moments), row.names = FALSE)
  cat("\n")

  cat("Dry and wet day counts:\n")
  print(.make_drywet_table(x$drydays), row.names = FALSE)
  cat("\n")

  cat("Selected percentiles (wet days only):\n")
  print(.make_quantiles_table(x$quantiles), row.names = FALSE)
  cat("\n")

  # --------------------------------------------------------------------------
  # Extremes (safe + label from actual thresholds)
  # --------------------------------------------------------------------------
  cat("Extreme precipitation changes (distributional comparison):\n")

  ex <- x$extremes
  if (is.null(ex) || !is.data.frame(ex) || nrow(ex) == 0) {
    cat("  (no extremes table available)\n\n")
  } else {

    if (!("threshold" %in% names(ex))) {
      cat("  (extremes table missing 'threshold' column)\n\n")
    } else {

      # Build readable labels from threshold values (e.g., 0.90 -> Top 10% (Percentile 90))
      # Interpret threshold as percentile probability p (0-1).
      p <- ex$threshold
      pct <- round(100 * p, 1)
      top <- round(100 * (1 - p), 3)

      # Pretty printing for common values
      fmt_pctile <- function(v) {
        ifelse(abs(v - round(v)) < 1e-8, sprintf("%.0f", round(v)), sprintf("%.1f", v))
      }
      fmt_top <- function(v) {
        # show 0.1% not 0.100%
        ifelse(v < 1, sprintf("%.1f%%", v), sprintf("%.0f%%", v))
      }

      ex$exceedance_label <- paste0(
        "Top ", fmt_top(top), " (Percentile ", fmt_pctile(pct), ")"
      )

      ext_print <- data.frame(
        `Exceedance level` = ex$exceedance_label,
        `Extreme cutoff change` = .format_pct(100 * (ex$threshold_ratio - 1), 0),
        `Average extreme size change` = .format_pct(100 * (ex$mean_ratio - 1), 0),
        `Maximum event change` = .format_pct(100 * (ex$max_ratio - 1), 0),
        stringsAsFactors = FALSE
      )

      print(ext_print, row.names = FALSE)
      cat("\n")
    }
  }

  cat("Spell lengths:\n")
  print(.make_spells_table(x$spells), row.names = FALSE)
  cat("\n")

  mtab <- .make_monthly_table(x$monthly)
  if (!is.null(mtab)) {
    cat("Monthly mean change vs target:\n")
    print(mtab, row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}






#' Summarize precipitation QM diagnostics
#'
#' @description
#' Prints a compact, numeric summary (no qualitative categories).
#'
#' @param object A `prcp_qm_diagnostics` object.
#' @param ... Unused.
#'
#' @export
summary.prcp_qm_diagnostics <- function(object, ...) {

  cat("\n=== Precipitation QM Diagnostics Summary ===\n\n")

  s <- object$summary

  cat("Moments:\n")
  cat(sprintf("  Mean ratio observed:   %.3f\n", s$moments$mean_ratio_observed))
  cat(sprintf("  Mean ratio intended:   %s\n", ifelse(is.finite(s$moments$mean_ratio_intended), sprintf("%.3f", s$moments$mean_ratio_intended), "NA")))
  cat(sprintf("  Mean ratio error:      %s\n", ifelse(is.finite(s$moments$mean_ratio_error), sprintf("%.3f", s$moments$mean_ratio_error), "NA")))
  cat(sprintf("  Mean within tolerance: %s\n", ifelse(is.na(s$moments$mean_within_tolerance), "NA", as.character(s$moments$mean_within_tolerance))))
  cat("\n")

  cat(sprintf("  Variance ratio observed:   %.3f\n", s$moments$var_ratio_observed))
  cat(sprintf("  Variance ratio intended:   %s\n", ifelse(is.finite(s$moments$var_ratio_intended), sprintf("%.3f", s$moments$var_ratio_intended), "NA")))
  cat(sprintf("  Variance ratio error:      %s\n", ifelse(is.finite(s$moments$var_ratio_error), sprintf("%.3f", s$moments$var_ratio_error), "NA")))
  cat(sprintf("  Variance within tolerance: %s\n", ifelse(is.na(s$moments$var_within_tolerance), "NA", as.character(s$moments$var_within_tolerance))))
  cat("\n")

  cat("Dry/Wet days (intended preserved):\n")
  cat(sprintf("  dry_days diff: %s; within_tolerance: %s\n",
              as.character(s$drydays$dry_days_diff),
              as.character(s$drydays$dry_days_within_tolerance)))
  cat(sprintf("  wet_days diff: %s; within_tolerance: %s\n",
              as.character(s$drydays$wet_days_diff),
              as.character(s$drydays$wet_days_within_tolerance)))
  cat("\n")

  cat("Quantiles (change magnitude):\n")
  cat(sprintf("  max |pct_change| across probs: %.2f%%\n", s$quantiles$max_abs_pct_change))
  cat("\n")

  cat("Extremes (ratios):\n")
  cat(sprintf("  mean threshold_ratio: %.3f\n", s$extremes$mean_threshold_ratio))
  cat(sprintf("  mean tail_mean_ratio: %.3f\n", s$extremes$mean_tail_mean_ratio))
  cat(sprintf("  mean max_ratio:       %.3f\n", s$extremes$mean_max_ratio))
  cat("\n")

  invisible(object)
}



#' Summarize precipitation QM diagnostics
#'
#' @description
#' Prints a compact numeric summary (no qualitative categories).
#'
#' @param object A `prcp_qm_diagnostics` object.
#' @param ... Unused.
#'
#' @export
summary.prcp_qm_diagnostics <- function(object, ...) {

  cat("\n=== Precipitation QM Diagnostics Summary ===\n\n")

  s <- object$summary

  # Moments (scenario-aware)
  cat("Overall targets vs observed:\n")
  if (is.finite(s$moments$mean_ratio_observed)) {
    cat(sprintf("  Mean:     observed %.3f", s$moments$mean_ratio_observed))
    if (is.finite(s$moments$mean_ratio_intended)) {
      cat(sprintf(" | intended %.3f | error %+0.3f\n",
                  s$moments$mean_ratio_intended,
                  s$moments$mean_ratio_error))
    } else {
      cat("\n")
    }
  }

  if (is.finite(s$moments$var_ratio_observed)) {
    cat(sprintf("  Variance: observed %.3f", s$moments$var_ratio_observed))
    if (is.finite(s$moments$var_ratio_intended)) {
      cat(sprintf(" | intended %.3f | error %+0.3f\n",
                  s$moments$var_ratio_intended,
                  s$moments$var_ratio_error))
    } else {
      cat("\n")
    }
  }
  cat("\n")

  # Dry/wet days (intended preserve)
  cat("Dry/wet day preservation (intended: no change):\n")
  cat(sprintf("  Dry days change: %s\n", as.character(s$drydays$dry_days_diff)))
  cat(sprintf("  Wet days change: %s\n", as.character(s$drydays$wet_days_diff)))
  cat("\n")

  # Quantiles (magnitude)
  cat("Quantiles (change magnitude):\n")
  cat(sprintf("  Max |pct change| across reported percentiles: %.2f%%\n",
              s$quantiles$max_abs_pct_change))
  cat("\n")

  # Extremes (distributional ratios)
  cat("Extremes (distributional ratios):\n")
  cat(sprintf("  Mean cutoff ratio:     %.3f\n", s$extremes$mean_threshold_ratio))
  cat(sprintf("  Mean tail-mean ratio:  %.3f\n", s$extremes$mean_tail_mean_ratio))
  cat(sprintf("  Mean max ratio:        %.3f\n", s$extremes$mean_max_ratio))
  cat("\n")

  invisible(object)
}


