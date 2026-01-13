# ==============================================================================
# VALIDATION DIAGNOSTICS FOR QUANTILE MAPPING
# ==============================================================================

#' Compute Validation Diagnostics for Quantile Mapping
#'
#' @description
#' Calculates comprehensive validation metrics comparing original and adjusted
#' precipitation time series. Evaluates moment preservation, distributional
#' characteristics, extreme values, and temporal patterns.
#'
#' @param value.original Numeric vector. Original precipitation values.
#' @param value.adjusted Numeric vector. Adjusted precipitation values (same length).
#' @param mon.ts Integer vector. Month indices (1-12) for each value.
#' @param year.ts Integer vector. Year indices for each value.
#' @param mean.change Numeric matrix. Target mean change factors (n_years x 12).
#' @param var.change Numeric matrix. Target variance change factors (n_years x 12).
#' @param wet.threshold Numeric. Threshold for defining wet days (default = 0.1 mm).
#' @param quantiles Numeric vector. Quantiles to evaluate (default = standard set).
#'
#' @return A list of class "qmap_diagnostics" containing:
#' \itemize{
#'   \item \code{moments}: Data frame of moment statistics
#'   \item \code{quantiles}: Data frame of quantile comparisons
#'   \item \code{extremes}: Data frame of extreme value metrics
#'   \item \code{temporal}: Data frame of temporal pattern metrics
#'   \item \code{monthly}: List of monthly-specific diagnostics
#'   \item \code{spells}: Data frame of spell length statistics
#'   \item \code{summary}: Overall assessment metrics
#' }
#'
#' @importFrom gridExtra grid.arrange
#' @export
validate_quantile_mapping <- function(
    value.original,
    value.adjusted,
    mon.ts = NULL,
    year.ts = NULL,
    mean.change = NULL,
    var.change = NULL,
    wet.threshold = 0.1,
    quantiles = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)) {

  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================

  if (length(value.original) != length(value.adjusted)) {
    stop("'value.original' and 'value.adjusted' must have same length")
  }

  n <- length(value.original)

  if (!is.null(mon.ts) && length(mon.ts) != n) {
    stop("'mon.ts' must have same length as values")
  }

  if (!is.null(year.ts) && length(year.ts) != n) {
    stop("'year.ts' must have same length as values")
  }

  # ==========================================================================
  # PREPARE DATA
  # ==========================================================================

  # Create masks for different subsets
  nonzero.original <- value.original > 0
  nonzero.adjusted <- value.adjusted > 0
  nonzero.both <- nonzero.original & nonzero.adjusted

  wet.original <- value.original >= wet.threshold
  wet.adjusted <- value.adjusted >= wet.threshold

  # ==========================================================================
  # 1. MOMENT DIAGNOSTICS
  # ==========================================================================

  moments <- compute_moment_diagnostics(
    value.original = value.original,
    value.adjusted = value.adjusted,
    nonzero.mask = nonzero.both
  )

  # ==========================================================================
  # 2. QUANTILE DIAGNOSTICS
  # ==========================================================================

  quantile.metrics <- compute_quantile_diagnostics(
    value.original = value.original[nonzero.both],
    value.adjusted = value.adjusted[nonzero.both],
    quantiles = quantiles
  )

  # ==========================================================================
  # 3. EXTREME VALUE DIAGNOSTICS
  # ==========================================================================

  extreme.metrics <- compute_extreme_diagnostics(
    value.original = value.original[nonzero.both],
    value.adjusted = value.adjusted[nonzero.both]
  )

  # ==========================================================================
  # 4. TEMPORAL PATTERN DIAGNOSTICS
  # ==========================================================================

  temporal.metrics <- NULL
  if (!is.null(mon.ts) && !is.null(year.ts)) {
    temporal.metrics <- compute_temporal_diagnostics(
      value.original = value.original,
      value.adjusted = value.adjusted,
      mon.ts = mon.ts,
      year.ts = year.ts,
      nonzero.mask = nonzero.both
    )
  }

  # ==========================================================================
  # 5. MONTHLY DIAGNOSTICS
  # ==========================================================================

  monthly.metrics <- NULL
  if (!is.null(mon.ts)) {
    monthly.metrics <- compute_monthly_diagnostics(
      value.original = value.original,
      value.adjusted = value.adjusted,
      mon.ts = mon.ts,
      mean.change = mean.change,
      var.change = var.change,
      year.ts = year.ts
    )
  }

  # ==========================================================================
  # 6. SPELL LENGTH DIAGNOSTICS
  # ==========================================================================

  spell.metrics <- compute_spell_diagnostics(
    value.original = value.original,
    value.adjusted = value.adjusted,
    threshold = wet.threshold
  )

  # ==========================================================================
  # 7. DRY DAY DIAGNOSTICS
  # ==========================================================================

  dryday.metrics <- compute_dryday_diagnostics(
    value.original = value.original,
    value.adjusted = value.adjusted,
    threshold = wet.threshold
  )

  # ==========================================================================
  # 8. OVERALL SUMMARY
  # ==========================================================================

  summary.metrics <- compute_summary_metrics(
    moments = moments,
    quantile.metrics = quantile.metrics,
    extreme.metrics = extreme.metrics,
    spell.metrics = spell.metrics,
    dryday.metrics = dryday.metrics
  )

  # ==========================================================================
  # RETURN RESULTS
  # ==========================================================================

  result <- list(
    moments = moments,
    quantiles = quantile.metrics,
    extremes = extreme.metrics,
    temporal = temporal.metrics,
    monthly = monthly.metrics,
    spells = spell.metrics,
    drydays = dryday.metrics,
    summary = summary.metrics
  )

  class(result) <- c("qmap_diagnostics", "list")
  result
}


# ==============================================================================
# HELPER FUNCTIONS FOR DIAGNOSTICS
# ==============================================================================

#' Compute moment diagnostics for quantile-mapped precipitation
#'
#' @description
#' Computes distribution moment statistics for original and adjusted series over a
#' user-specified subset (typically days where both series are non-zero).
#' Returns absolute and relative changes and assigns a qualitative assessment per metric.
#'
#' @param value.original Numeric vector. Original precipitation values.
#' @param value.adjusted Numeric vector. Adjusted precipitation values. Must be the same length
#'   as \code{value.original}.
#' @param nonzero.mask Logical vector. Subset mask applied to both series before computing
#'   statistics (e.g., \code{value.original > 0 & value.adjusted > 0}).
#'
#' @details
#' Metrics computed on the masked series are: mean, standard deviation, variance,
#' coefficient of variation (CV), skewness, and excess kurtosis. Ratio, difference,
#' and percent change are computed as \code{adjusted/original}, \code{adjusted-original},
#' and \code{(ratio - 1) * 100}.
#'
#' @return A data.frame with one row per metric and columns:
#' \itemize{
#'   \item \code{metric}: Metric name.
#'   \item \code{original}, \code{adjusted}: Metric values.
#'   \item \code{ratio}, \code{diff}, \code{pct.change}: Relative/absolute change diagnostics.
#'   \item \code{assessment}: Qualitative assessment (see \code{\link{assess_moment_changes}}).
#' }
#'
#' @keywords internal
compute_moment_diagnostics <- function(value.original, value.adjusted, nonzero.mask) {

  # Overall statistics
  overall <- data.frame(
    metric = c("mean", "sd", "variance", "cv", "skewness", "kurtosis"),
    original = c(
      mean(value.original[nonzero.mask]),
      sd(value.original[nonzero.mask]),
      var(value.original[nonzero.mask]),
      sd(value.original[nonzero.mask]) / mean(value.original[nonzero.mask]),
      compute_skewness(value.original[nonzero.mask]),
      compute_kurtosis(value.original[nonzero.mask])
    ),
    adjusted = c(
      mean(value.adjusted[nonzero.mask]),
      sd(value.adjusted[nonzero.mask]),
      var(value.adjusted[nonzero.mask]),
      sd(value.adjusted[nonzero.mask]) / mean(value.adjusted[nonzero.mask]),
      compute_skewness(value.adjusted[nonzero.mask]),
      compute_kurtosis(value.adjusted[nonzero.mask])
    ),
    stringsAsFactors = FALSE
  )

  # Compute ratios and differences
  overall$ratio <- overall$adjusted / overall$original
  overall$diff <- overall$adjusted - overall$original
  overall$pct.change <- (overall$ratio - 1) * 100

  # Add assessment
  overall$assessment <- assess_moment_changes(overall)

  overall
}


#' Compute quantile diagnostics for quantile-mapped precipitation
#'
#' @description
#' Compares specified quantiles between original and adjusted series and reports
#' relative and absolute differences with a simple qualitative assessment.
#'
#' @param value.original Numeric vector. Original precipitation values (typically filtered
#'   to non-zero or wet days upstream).
#' @param value.adjusted Numeric vector. Adjusted precipitation values, same length as
#'   \code{value.original}.
#' @param quantiles Numeric vector in [0, 1]. Quantile probabilities to evaluate.
#'
#' @details
#' Quantiles are computed with \code{stats::quantile(..., na.rm = TRUE)}. Percent change
#' is computed from the ratio \code{adjusted/original}. Assessment thresholds currently are:
#' \itemize{
#'   \item \code{"good"}: |pct.change| < 5
#'   \item \code{"acceptable"}: 5 <= |pct.change| < 15
#'   \item \code{"poor"}: |pct.change| >= 15
#' }
#'
#' @return A data.frame with one row per requested probability and columns:
#' \code{quantile}, \code{original}, \code{adjusted}, \code{ratio}, \code{diff},
#' \code{pct.change}, \code{abs.error}, \code{assessment}.
#'
#' @keywords internal
compute_quantile_diagnostics <- function(value.original, value.adjusted, quantiles) {

  q.original <- quantile(value.original, probs = quantiles, na.rm = TRUE)
  q.adjusted <- quantile(value.adjusted, probs = quantiles, na.rm = TRUE)

  result <- data.frame(
    quantile = quantiles,
    original = as.numeric(q.original),
    adjusted = as.numeric(q.adjusted),
    stringsAsFactors = FALSE
  )

  result$ratio <- result$adjusted / result$original
  result$diff <- result$adjusted - result$original
  result$pct.change <- (result$ratio - 1) * 100

  # Compute RMSE and MAE
  result$abs.error <- abs(result$diff)

  # Add assessment
  result$assessment <- ifelse(
    abs(result$pct.change) < 5, "good",
    ifelse(abs(result$pct.change) < 15, "acceptable", "poor")
  )

  result
}


#' Compute extreme-value diagnostics for quantile-mapped precipitation
#'
#' @description
#' Evaluates tail behavior by comparing exceedance thresholds and summary statistics
#' of values above those thresholds for original and adjusted series.
#'
#' @param value.original Numeric vector. Original precipitation values (typically filtered
#'   to non-zero or wet days upstream).
#' @param value.adjusted Numeric vector. Adjusted precipitation values, same length as
#'   \code{value.original}.
#'
#' @details
#' For each threshold probability \code{p} in \code{c(0.90, 0.95, 0.99, 0.999)},
#' the function computes the empirical quantile \code{Q(p)} and then computes
#' the mean and maximum of the subset \code{x[x >= Q(p)]}. Ratios are returned for
#' threshold, mean, and max.
#'
#' @return A data.frame with one row per threshold probability and columns:
#' \itemize{
#'   \item \code{threshold}: Probability used for \code{quantile()}.
#'   \item \code{original.threshold}, \code{adjusted.threshold}: Exceedance cutoffs.
#'   \item \code{original.mean}, \code{adjusted.mean}: Mean above cutoff.
#'   \item \code{original.max}, \code{adjusted.max}: Max above cutoff.
#'   \item \code{threshold.ratio}, \code{mean.ratio}, \code{max.ratio}: Adjusted/original ratios.
#' }
#'
#' @keywords internal
compute_extreme_diagnostics <- function(value.original, value.adjusted) {

  # Define extreme thresholds
  thresholds <- c(0.90, 0.95, 0.99, 0.999)

  metrics <- data.frame(
    threshold = thresholds,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]

    # Original extremes
    q.orig <- quantile(value.original, thresh, na.rm = TRUE)
    extremes.orig <- value.original[value.original >= q.orig]

    # Adjusted extremes
    q.adj <- quantile(value.adjusted, thresh, na.rm = TRUE)
    extremes.adj <- value.adjusted[value.adjusted >= q.adj]

    # Statistics
    metrics$original.threshold[i] <- q.orig
    metrics$adjusted.threshold[i] <- q.adj
    metrics$original.mean[i] <- mean(extremes.orig, na.rm = TRUE)
    metrics$adjusted.mean[i] <- mean(extremes.adj, na.rm = TRUE)
    metrics$original.max[i] <- max(extremes.orig, na.rm = TRUE)
    metrics$adjusted.max[i] <- max(extremes.adj, na.rm = TRUE)
  }

  # Compute ratios
  metrics$threshold.ratio <- metrics$adjusted.threshold / metrics$original.threshold
  metrics$mean.ratio <- metrics$adjusted.mean / metrics$original.mean
  metrics$max.ratio <- metrics$adjusted.max / metrics$original.max

  metrics
}


#' Compute temporal diagnostics for quantile-mapped precipitation
#'
#' @description
#' Quantifies preservation of temporal structure by comparing annual and monthly means,
#' their correlations, and lag-1 autocorrelation of the (masked) daily series.
#'
#' @param value.original Numeric vector. Original precipitation values.
#' @param value.adjusted Numeric vector. Adjusted precipitation values, same length as
#'   \code{value.original}.
#' @param mon.ts Integer vector. Month index (1--12) for each observation.
#' @param year.ts Integer vector. Year index for each observation (same length as values).
#' @param nonzero.mask Logical vector. Mask applied to both series before aggregation
#'   and autocorrelation calculations.
#'
#' @details
#' Annual and monthly means are computed via \code{tapply(..., mean, na.rm = TRUE)} on
#' the masked values. Correlations compare original vs adjusted aggregated series.
#' Lag-1 autocorrelation is computed using \code{stats::acf(..., lag.max = 1, plot = FALSE)}.
#'
#' @return A data.frame with columns \code{metric}, \code{value}, and \code{assessment}.
#' The \code{assessment} field is populated for correlation and autocorrelation-preservation
#' diagnostics and set to \code{NA} for raw autocorrelation values.
#'
#' @keywords internal
compute_temporal_diagnostics <- function(value.original, value.adjusted,
                                         mon.ts, year.ts, nonzero.mask) {

  # Annual means
  annual.orig <- tapply(
    value.original[nonzero.mask],
    year.ts[nonzero.mask],
    mean, na.rm = TRUE
  )

  annual.adj <- tapply(
    value.adjusted[nonzero.mask],
    year.ts[nonzero.mask],
    mean, na.rm = TRUE
  )

  # Monthly means
  monthly.orig <- tapply(
    value.original[nonzero.mask],
    mon.ts[nonzero.mask],
    mean, na.rm = TRUE
  )

  monthly.adj <- tapply(
    value.adjusted[nonzero.mask],
    mon.ts[nonzero.mask],
    mean, na.rm = TRUE
  )

  # Compute correlations
  cor.annual <- cor(annual.orig, annual.adj, use = "complete.obs")
  cor.monthly <- cor(monthly.orig, monthly.adj, use = "complete.obs")

  # Compute autocorrelations
  acf.orig <- acf(value.original[nonzero.mask], lag.max = 1, plot = FALSE)$acf[2]
  acf.adj <- acf(value.adjusted[nonzero.mask], lag.max = 1, plot = FALSE)$acf[2]

  data.frame(
    metric = c(
      "annual.correlation",
      "monthly.correlation",
      "lag1.autocorr.original",
      "lag1.autocorr.adjusted",
      "autocorr.preservation"
    ),
    value = c(
      cor.annual,
      cor.monthly,
      acf.orig,
      acf.adj,
      abs(acf.adj - acf.orig)
    ),
    assessment = c(
      ifelse(cor.annual > 0.9, "excellent", ifelse(cor.annual > 0.7, "good", "poor")),
      ifelse(cor.monthly > 0.95, "excellent", ifelse(cor.monthly > 0.85, "good", "poor")),
      NA, NA,
      ifelse(abs(acf.adj - acf.orig) < 0.1, "good", "moderate")
    ),
    stringsAsFactors = FALSE
  )
}


#' Compute month-by-month diagnostics for quantile-mapped precipitation
#'
#' @description
#' Computes monthly mean/variance diagnostics for original and adjusted precipitation,
#' optionally comparing achieved changes against target mean/variance change factors.
#'
#' @param value.original Numeric vector. Original precipitation values.
#' @param value.adjusted Numeric vector. Adjusted precipitation values, same length as
#'   \code{value.original}.
#' @param mon.ts Integer vector. Month index (1--12) for each observation.
#' @param mean.change Optional numeric matrix (n_years x 12). Target mean change factors
#'   by year and month. If provided with \code{var.change} and \code{year.ts}, the function
#'   computes mean absolute error to targets.
#' @param var.change Optional numeric matrix (n_years x 12). Target variance change factors
#'   by year and month.
#' @param year.ts Optional integer vector. Year index for each observation. Required only
#'   when comparing to \code{mean.change} and \code{var.change}.
#'
#' @details
#' Monthly statistics are computed on indices where both series are strictly positive
#' (\code{value.original > 0 & value.adjusted > 0}) for each month. This isolates intensity
#' behavior from frequency effects. If target change matrices are provided, targets are
#' averaged across years using \code{colMeans(..., na.rm = TRUE)}.
#'
#' @return A data.frame with one row per month and columns including:
#' \code{original.mean}, \code{adjusted.mean}, \code{original.var}, \code{adjusted.var},
#' \code{n.obs}, \code{mean.ratio}, \code{var.ratio}. If targets are provided,
#' additional columns \code{target.mean.ratio}, \code{target.var.ratio}, \code{mean.error},
#' \code{var.error}, \code{mean.assessment}, \code{var.assessment} are included.
#'
#' @keywords internal
compute_monthly_diagnostics <- function(value.original, value.adjusted, mon.ts,
                                        mean.change = NULL, var.change = NULL,
                                        year.ts = NULL) {

  monthly <- data.frame(
    month = 1:12,
    stringsAsFactors = FALSE
  )

  for (m in 1:12) {
    idx <- mon.ts == m & value.original > 0 & value.adjusted > 0

    if (sum(idx) > 0) {
      monthly$original.mean[m] <- mean(value.original[idx], na.rm = TRUE)
      monthly$adjusted.mean[m] <- mean(value.adjusted[idx], na.rm = TRUE)
      monthly$original.var[m] <- var(value.original[idx], na.rm = TRUE)
      monthly$adjusted.var[m] <- var(value.adjusted[idx], na.rm = TRUE)
      monthly$n.obs[m] <- sum(idx)
    } else {
      monthly$original.mean[m] <- NA
      monthly$adjusted.mean[m] <- NA
      monthly$original.var[m] <- NA
      monthly$adjusted.var[m] <- NA
      monthly$n.obs[m] <- 0
    }
  }

  monthly$mean.ratio <- monthly$adjusted.mean / monthly$original.mean
  monthly$var.ratio <- monthly$adjusted.var / monthly$original.var

  # Compare with target changes if provided
  if (!is.null(mean.change) && !is.null(var.change) && !is.null(year.ts)) {
    # Average target changes across years
    monthly$target.mean.ratio <- colMeans(mean.change, na.rm = TRUE)
    monthly$target.var.ratio <- colMeans(var.change, na.rm = TRUE)

    monthly$mean.error <- abs(monthly$mean.ratio - monthly$target.mean.ratio)
    monthly$var.error <- abs(monthly$var.ratio - monthly$target.var.ratio)

    monthly$mean.assessment <- ifelse(
      monthly$mean.error < 0.05, "excellent",
      ifelse(monthly$mean.error < 0.15, "good", "poor")
    )

    monthly$var.assessment <- ifelse(
      monthly$var.error < 0.10, "excellent",
      ifelse(monthly$var.error < 0.25, "good", "poor")
    )
  }

  monthly
}


#' Compute wet and dry spell-length diagnostics
#'
#' @description
#' Summarizes persistence characteristics by comparing dry- and wet-spell lengths between
#' original and adjusted series using a user-defined wet-day threshold.
#'
#' @param value.original Numeric vector. Original precipitation values.
#' @param value.adjusted Numeric vector. Adjusted precipitation values, same length as
#'   \code{value.original}.
#' @param threshold Numeric scalar. Wet-day threshold; values < threshold are classified as dry
#'   and values >= threshold as wet.
#'
#' @details
#' Spell lengths are computed using \code{\link{compute_spell_lengths}} for dry spells
#' (\code{below = TRUE}) and wet spells (\code{below = FALSE}). The function reports mean,
#' maximum, and standard deviation of spell lengths and their ratios (adjusted/original).
#'
#' @return A data.frame with two rows (\code{"dry"}, \code{"wet"}) and columns:
#' \code{original.mean}, \code{adjusted.mean}, \code{original.max}, \code{adjusted.max},
#' \code{original.sd}, \code{adjusted.sd}, \code{mean.ratio}, \code{max.ratio}, \code{sd.ratio},
#' and \code{assessment}.
#'
#' @keywords internal
compute_spell_diagnostics <- function(value.original, value.adjusted, threshold) {

  # Dry spells
  dry.orig <- compute_spell_lengths(value.original, threshold, below = TRUE)
  dry.adj <- compute_spell_lengths(value.adjusted, threshold, below = TRUE)

  # Wet spells
  wet.orig <- compute_spell_lengths(value.original, threshold, below = FALSE)
  wet.adj <- compute_spell_lengths(value.adjusted, threshold, below = FALSE)

  result <- data.frame(
    spell.type = c("dry", "wet"),
    original.mean = c(mean(dry.orig), mean(wet.orig)),
    adjusted.mean = c(mean(dry.adj), mean(wet.adj)),
    original.max = c(max(dry.orig), max(wet.orig)),
    adjusted.max = c(max(dry.adj), max(wet.adj)),
    original.sd = c(sd(dry.orig), sd(wet.orig)),
    adjusted.sd = c(sd(dry.adj), sd(wet.adj)),
    stringsAsFactors = FALSE
  )

  result$mean.ratio <- result$adjusted.mean / result$original.mean
  result$max.ratio <- result$adjusted.max / result$original.max
  result$sd.ratio <- result$adjusted.sd / result$original.sd

  # Assessment: spell lengths should be relatively preserved
  result$assessment <- ifelse(
    abs(result$mean.ratio - 1) < 0.2, "good",
    ifelse(abs(result$mean.ratio - 1) < 0.4, "acceptable", "poor")
  )

  result
}


#' Compute dry-day and wet-day frequency diagnostics
#'
#' @description
#' Compares dry and wet day counts and frequencies between original and adjusted series
#' given a wet-day threshold.
#'
#' @param value.original Numeric vector. Original precipitation values.
#' @param value.adjusted Numeric vector. Adjusted precipitation values, same length as
#'   \code{value.original}.
#' @param threshold Numeric scalar. Wet-day threshold; values < threshold are classified as dry
#'   and values >= threshold as wet.
#'
#' @details
#' The function reports both counts and frequencies (count divided by total length).
#' Percent change is computed as \code{adjusted/original - 1}. If original counts are zero
#' (edge case), percent change may be \code{Inf} or \code{NaN}.
#'
#' @return A data.frame with four rows (dry/wet counts and frequencies) and columns:
#' \code{category}, \code{original}, \code{adjusted}, \code{difference}, \code{pct.change},
#' \code{assessment}.
#'
#' @keywords internal
compute_dryday_diagnostics <- function(value.original, value.adjusted, threshold) {

  n.total <- length(value.original)

  dry.orig <- sum(value.original < threshold)
  dry.adj <- sum(value.adjusted < threshold)

  wet.orig <- sum(value.original >= threshold)
  wet.adj <- sum(value.adjusted >= threshold)

  result <- data.frame(
    category = c("dry.days", "wet.days", "dry.frequency", "wet.frequency"),
    original = c(dry.orig, wet.orig, dry.orig / n.total, wet.orig / n.total),
    adjusted = c(dry.adj, wet.adj, dry.adj / n.total, wet.adj / n.total),
    stringsAsFactors = FALSE
  )

  result$difference <- result$adjusted - result$original
  result$pct.change <- (result$adjusted / result$original - 1) * 100

  # Assessment: dry day frequency should be exactly preserved
  result$assessment <- c(
    ifelse(result$difference[1] == 0, "perfect", "changed"),
    ifelse(result$difference[2] == 0, "perfect", "changed"),
    ifelse(abs(result$difference[3]) < 0.01, "excellent", "changed"),
    ifelse(abs(result$difference[4]) < 0.01, "excellent", "changed")
  )

  result
}


#' Compute overall summary metrics for quantile-mapping diagnostics
#'
#' @description
#' Aggregates component diagnostics into an overall score and qualitative assessment,
#' and extracts a compact set of key findings for reporting.
#'
#' @param moments Data.frame. Output from \code{\link{compute_moment_diagnostics}}.
#' @param quantile.metrics Data.frame. Output from \code{\link{compute_quantile_diagnostics}}.
#' @param extreme.metrics Data.frame. Output from \code{\link{compute_extreme_diagnostics}}.
#' @param spell.metrics Data.frame. Output from \code{\link{compute_spell_diagnostics}}.
#' @param dryday.metrics Data.frame. Output from \code{\link{compute_dryday_diagnostics}}.
#'
#' @details
#' Current scoring logic (0--100) weights:
#' moments (30), quantiles (30), extremes (20), spells (10), dry days (10).
#' The score is designed for quick triage and is not a formal statistical test.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{overall.score}: Numeric scalar in [0, 100] (after bounding where applied).
#'   \item \code{component.scores}: Named numeric vector with per-component scores.
#'   \item \code{overall.assessment}: Character scalar in \code{c("excellent","good","acceptable","poor")}.
#'   \item \code{key.findings}: List of headline diagnostics for reporting.
#' }
#'
#' @keywords internal
compute_summary_metrics <- function(moments, quantile.metrics, extreme.metrics,
                                    spell.metrics, dryday.metrics) {

  summary <- list()

  # Overall quality score (0-100)
  scores <- numeric()

  # Moment preservation (30 points)
  mean.score <- ifelse(abs(moments$pct.change[1]) < 5, 10,
                       ifelse(abs(moments$pct.change[1]) < 15, 5, 0))
  var.score <- ifelse(abs(moments$pct.change[3]) < 10, 10,
                      ifelse(abs(moments$pct.change[3]) < 25, 5, 0))
  cv.score <- ifelse(abs(moments$pct.change[4]) < 10, 10,
                     ifelse(abs(moments$pct.change[4]) < 20, 5, 0))
  scores["moments"] <- mean.score + var.score + cv.score

  # Quantile preservation (30 points)
  good.quantiles <- sum(quantile.metrics$assessment == "good")
  acceptable.quantiles <- sum(quantile.metrics$assessment == "acceptable")
  quantile.score <- (good.quantiles * 4 + acceptable.quantiles * 2)
  scores["quantiles"] <- min(quantile.score, 30)

  # Extreme value preservation (20 points)
  extreme.score <- 20 * (1 - mean(abs(extreme.metrics$mean.ratio - 1)))
  scores["extremes"] <- max(0, min(extreme.score, 20))

  # Spell preservation (10 points)
  spell.score <- sum(spell.metrics$assessment == "good") * 5
  scores["spells"] <- min(spell.score, 10)

  # Dry day preservation (10 points)
  dryday.score <- ifelse(dryday.metrics$assessment[1] == "perfect", 10, 0)
  scores["drydays"] <- dryday.score

  summary$overall.score <- sum(scores)
  summary$component.scores <- scores

  summary$overall.assessment <- ifelse(
    summary$overall.score >= 80, "excellent",
    ifelse(summary$overall.score >= 60, "good",
           ifelse(summary$overall.score >= 40, "acceptable", "poor"))
  )

  # Key findings
  summary$key.findings <- list(
    mean.change.pct = moments$pct.change[moments$metric == "mean"],
    var.change.pct = moments$pct.change[moments$metric == "variance"],
    dryday.preserved = dryday.metrics$difference[1] == 0,
    worst.quantile = quantile.metrics$quantile[which.max(abs(quantile.metrics$pct.change))],
    extreme.amplification = mean(extreme.metrics$mean.ratio)
  )

  summary
}

# ==============================================================================
# PRINT AND PLOT METHODS
# ==============================================================================

#' Print quantile-mapping validation diagnostics
#'
#' @description
#' S3 print method for objects of class \code{"qmap_diagnostics"}. Produces a structured,
#' human-readable console report summarizing overall performance, component scores,
#' key findings, and selected detailed diagnostics.
#'
#' @param x An object of class \code{"qmap_diagnostics"} as returned by
#'   \code{\link{validate_quantile_mapping}}.
#' @param ... Additional arguments (currently ignored). Included for S3 method compatibility.
#'
#' @details
#' The printed output includes:
#' \itemize{
#'   \item Overall assessment label and total score (0--100).
#'   \item Component-level scores for moments, quantiles, extremes, spells, and dry days.
#'   \item Key findings (mean/variance change, dry-day preservation, extreme amplification).
#'   \item Tabular summaries of moment preservation and dry-day counts.
#' }
#' This method is intended for interactive inspection. For programmatic access,
#' use the corresponding elements of the \code{"qmap_diagnostics"} object directly
#' (e.g., \code{x$summary}, \code{x$moments}, \code{x$drydays}).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @seealso \code{\link{validate_quantile_mapping}},
#'   \code{\link{summary.qmap_diagnostics}},
#'   \code{\link{plot.qmap_diagnostics}}
#'
#' @export
print.qmap_diagnostics <- function(x, ...) {
  cat("\n=== Quantile Mapping Diagnostics ===\n\n")

  # Overall summary
  cat("Overall Assessment:", x$summary$overall.assessment, "\n")
  cat("Overall Score:", round(x$summary$overall.score, 1), "/ 100\n\n")

  # Component scores
  cat("Component Scores:\n")
  for (comp in names(x$summary$component.scores)) {
    cat(sprintf("  %s: %.1f\n", comp, x$summary$component.scores[comp]))
  }
  cat("\n")

  # Key findings
  cat("Key Findings:\n")
  cat(sprintf("  Mean change: %.1f%%\n", x$summary$key.findings$mean.change.pct))
  cat(sprintf("  Variance change: %.1f%%\n", x$summary$key.findings$var.change.pct))
  cat(sprintf("  Dry days preserved: %s\n", x$summary$key.findings$dryday.preserved))
  cat(sprintf("  Extreme amplification: %.2fx\n", x$summary$key.findings$extreme.amplification))
  cat("\n")

  # Moment preservation
  cat("Moment Preservation:\n")
  print(x$moments[, c("metric", "original", "adjusted", "pct.change", "assessment")],
        row.names = FALSE, digits = 3)
  cat("\n")

  # Dry day preservation
  cat("Dry Day Preservation:\n")
  print(x$drydays[1:2, c("category", "original", "adjusted", "assessment")],
        row.names = FALSE)
  cat("\n")

  invisible(x)
}



#' Plot diagnostics for quantile-mapping validation
#'
#' @description
#' S3 plot method for objects of class \code{"qmap_diagnostics"} produced by
#' \code{\link{validate_quantile_mapping}}. Generates one or more ggplot-based diagnostic
#' figures summarizing moment changes, quantile preservation, extreme-value behavior,
#' and (optionally) monthly mean patterns.
#'
#' @param x An object of class \code{"qmap_diagnostics"} as returned by
#'   \code{\link{validate_quantile_mapping}}.
#' @param which Character scalar. Which plot to produce:
#' \itemize{
#'   \item \code{"all"}: return/arrange all available panels.
#'   \item \code{"moments"}: percent changes in distribution moments.
#'   \item \code{"quantiles"}: scatter of original vs adjusted quantiles.
#'   \item \code{"extremes"}: comparison of mean exceedances across thresholds.
#'   \item \code{"monthly"}: monthly mean comparison (requires \code{x$monthly}).
#' }
#' @param ... Additional arguments (currently ignored). Included for S3 method compatibility.
#'
#' @details
#' This method requires \pkg{ggplot2}. If \code{which = "monthly"}, it also requires
#' \pkg{tidyr} for reshaping. If \code{which = "all"} and \pkg{gridExtra} is installed,
#' plots are arranged in a multi-panel layout; otherwise a named list of ggplot objects
#' is returned.
#'
#' @return
#' If \code{which = "all"}:
#' \itemize{
#'   \item If \pkg{gridExtra} is available, the arranged grob is drawn and returned
#'   invisibly (side effect).
#'   \item Otherwise, a named list of ggplot objects is returned.
#' }
#' If \code{which != "all"}, a single ggplot object is returned.
#'
#' @seealso \code{\link{validate_quantile_mapping}}, \code{\link{print.qmap_diagnostics}},
#'   \code{\link{summary.qmap_diagnostics}}
#'
#' @export
plot.qmap_diagnostics <- function(x, which = c("all", "moments", "quantiles",
                                               "extremes", "monthly"), ...) {

  which <- match.arg(which)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }

  plots <- list()

  # 1. Moments comparison
  if (which %in% c("all", "moments")) {
    p1 <- ggplot2::ggplot(x$moments, ggplot2::aes(x = metric)) +
      ggplot2::geom_col(ggplot2::aes(y = pct.change, fill = assessment)) +
      ggplot2::scale_fill_manual(
        values = c("excellent" = "green3", "good" = "yellow3", "poor" = "red3")
      ) +
      ggplot2::labs(
        title = "Moment Changes",
        x = "Metric",
        y = "Percent Change (%)",
        fill = "Assessment"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    plots$moments <- p1
  }

  # 2. Quantile comparison
  if (which %in% c("all", "quantiles")) {
    p2 <- ggplot2::ggplot(x$quantiles, ggplot2::aes(x = original, y = adjusted)) +
      ggplot2::geom_point(ggplot2::aes(color = assessment), size = 3) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::scale_color_manual(
        values = c("good" = "green3", "acceptable" = "yellow3", "poor" = "red3")
      ) +
      ggplot2::labs(
        title = "Quantile Preservation",
        x = "Original Quantile Value",
        y = "Adjusted Quantile Value",
        color = "Assessment"
      ) +
      ggplot2::theme_minimal()

    plots$quantiles <- p2
  }

  # 3. Extreme values
  if (which %in% c("all", "extremes")) {
    extreme.long <- data.frame(
      threshold = rep(x$extremes$threshold, 2),
      type = rep(c("Original", "Adjusted"), each = nrow(x$extremes)),
      mean = c(x$extremes$original.mean, x$extremes$adjusted.mean)
    )

    p3 <- ggplot2::ggplot(extreme.long, ggplot2::aes(
      x = factor(threshold),
      y = mean,
      fill = type
    )) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::scale_fill_manual(values = c("Original" = "blue3", "Adjusted" = "red3")) +
      ggplot2::labs(
        title = "Extreme Value Comparison",
        x = "Exceedance Quantile",
        y = "Mean Value",
        fill = ""
      ) +
      ggplot2::theme_minimal()

    plots$extremes <- p3
  }

  # 4. Monthly patterns
  if (which %in% c("all", "monthly") && !is.null(x$monthly)) {
    monthly.long <- tidyr::pivot_longer(
      x$monthly[, c("month", "original.mean", "adjusted.mean")],
      cols = c("original.mean", "adjusted.mean"),
      names_to = "type",
      values_to = "mean"
    )
    monthly.long$type <- ifelse(
      monthly.long$type == "original.mean", "Original", "Adjusted"
    )

    p4 <- ggplot2::ggplot(monthly.long, ggplot2::aes(
      x = factor(month),
      y = mean,
      color = type,
      group = type
    )) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = c("Original" = "blue3", "Adjusted" = "red3")) +
      ggplot2::scale_x_discrete(labels = month.abb) +
      ggplot2::labs(
        title = "Monthly Mean Precipitation",
        x = "Month",
        y = "Mean (mm/day)",
        color = ""
      ) +
      ggplot2::theme_minimal()

    plots$monthly <- p4
  }

  # Return plots
  if (which == "all") {
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(grobs = plots, ncol = 2)
    } else {
      return(plots)
    }
  } else {
    return(plots[[1]])
  }
}


#' Summarize quantile-mapping validation diagnostics
#'
#' @description
#' S3 summary method for \code{"qmap_diagnostics"} objects. Prints a concise textual
#' summary including the overall score/assessment and selected component-level
#' diagnostics (moments, quantiles, dry days, and extremes).
#'
#' @param object An object of class \code{"qmap_diagnostics"} as returned by
#'   \code{\link{validate_quantile_mapping}}.
#' @param ... Additional arguments (currently ignored). Included for S3 method compatibility.
#'
#' @details
#' This method writes to the console via \code{cat()}. It is intended for quick
#' human-readable reporting rather than programmatic extraction (use \code{object$summary}
#' and related list elements for that).
#'
#' @return The input object, invisibly.
#'
#' @seealso \code{\link{validate_quantile_mapping}}, \code{\link{print.qmap_diagnostics}},
#'   \code{\link{plot.qmap_diagnostics}}
#'
#' @export
summary.qmap_diagnostics <- function(object, ...) {
  cat("\n=== Quantile Mapping Validation Summary ===\n\n")

  cat("Overall Score:", round(object$summary$overall.score, 1), "/ 100\n")
  cat("Assessment:", object$summary$overall.assessment, "\n\n")

  cat("Detailed Assessment:\n")
  cat("------------------\n")

  # Moments
  cat("\nMoments:\n")
  good.moments <- sum(object$moments$assessment %in% c("excellent", "good"))
  cat(sprintf("  %d/%d metrics within acceptable range\n",
              good.moments, nrow(object$moments)))

  # Quantiles
  cat("\nQuantiles:\n")
  good.quantiles <- sum(object$quantiles$assessment == "good")
  cat(sprintf("  %d/%d quantiles well preserved\n",
              good.quantiles, nrow(object$quantiles)))

  # Dry days
  cat("\nDry Days:\n")
  cat(sprintf("  Preservation: %s\n", object$drydays$assessment[1]))

  # Extremes
  cat("\nExtremes:\n")
  cat(sprintf("  Mean amplification: %.2fx\n",
              object$summary$key.findings$extreme.amplification))

  invisible(object)
}
