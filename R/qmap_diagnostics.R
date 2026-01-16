# ==============================================================================
# VALIDATION DIAGNOSTICS FOR QUANTILE MAPPING
# ==============================================================================

#' Compute Diagnostics for Precipitation Quantile Mapping
#'
#' @description
#' Calculates validation-style diagnostics comparing reference and quantile-mapped
#' precipitation series. Evaluates moment changes, quantile preservation, tail behavior,
#' dry/wet-day frequencies, spell lengths, and optional temporal/monthly patterns.
#'
#' @param prcp_ref Numeric vector. Reference (original) precipitation values.
#' @param prcp_adj Numeric vector. Adjusted precipitation values (same length as \code{prcp_ref}).
#' @param month Integer vector (optional). Month index (1--12) for each observation.
#' @param year Integer vector (optional). Year index (1..n_years) for each observation.
#' @param mean_factor Numeric matrix (optional). Target mean scaling factors (n_years x 12).
#' @param var_factor Numeric matrix (optional). Target variance scaling factors (n_years x 12).
#' @param wet_thresh Numeric. Threshold for defining wet days (default = 0.1 mm).
#' @param probs Numeric vector. Quantile probabilities to evaluate.
#'
#' @return A list of class \code{"prcp_qm_diagnostics"} containing:
#' \itemize{
#'   \item \code{moments}: Data frame of moment statistics
#'   \item \code{quantiles}: Data frame of quantile comparisons
#'   \item \code{extremes}: Data frame of extreme value metrics
#'   \item \code{temporal}: Data frame of temporal pattern metrics (optional)
#'   \item \code{monthly}: Data frame of monthly diagnostics (optional)
#'   \item \code{spells}: Data frame of spell length statistics
#'   \item \code{drydays}: Data frame of dry/wet-day frequency diagnostics
#'   \item \code{summary}: Overall assessment metrics
#' }
#'
#' @importFrom gridExtra grid.arrange
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

  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================

  if (length(prcp_ref) != length(prcp_adj)) {
    stop("'prcp_ref' and 'prcp_adj' must have same length")
  }

  n <- length(prcp_ref)

  if (!is.null(month) && length(month) != n) {
    stop("'month' must have same length as precipitation series")
  }

  if (!is.null(year) && length(year) != n) {
    stop("'year' must have same length as precipitation series")
  }

  # ==========================================================================
  # PREPARE DATA
  # ==========================================================================

  # Masks
  is_nz_ref <- prcp_ref > 0
  is_nz_adj <- prcp_adj > 0
  is_nz_both <- is_nz_ref & is_nz_adj

  is_wetday_ref <- prcp_ref >= wet_thresh
  is_wetday_adj <- prcp_adj >= wet_thresh
  is_wetday_both <- is_wetday_ref & is_wetday_adj

  # ==========================================================================
  # 1. MOMENT DIAGNOSTICS
  # ==========================================================================

  moments <- compute_moment_diagnostics(
    prcp_ref = prcp_ref,
    prcp_adj = prcp_adj,
    mask = is_nz_both
  )

  # ==========================================================================
  # 2. QUANTILE DIAGNOSTICS
  # ==========================================================================

  quant_metrics <- compute_quantile_diagnostics(
    prcp_ref = prcp_ref[is_nz_both],
    prcp_adj = prcp_adj[is_nz_both],
    probs = probs
  )

  # ==========================================================================
  # 3. EXTREME VALUE DIAGNOSTICS
  # ==========================================================================

  extreme_metrics <- compute_extreme_diagnostics(
    prcp_ref = prcp_ref[is_nz_both],
    prcp_adj = prcp_adj[is_nz_both]
  )

  # ==========================================================================
  # 4. TEMPORAL PATTERN DIAGNOSTICS
  # ==========================================================================

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

  # ==========================================================================
  # 5. MONTHLY DIAGNOSTICS
  # ==========================================================================

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

  # ==========================================================================
  # 6. SPELL LENGTH DIAGNOSTICS
  # ==========================================================================

  spell_metrics <- compute_spell_diagnostics(
    prcp_ref = prcp_ref,
    prcp_adj = prcp_adj,
    wet_thresh = wet_thresh
  )

  # ==========================================================================
  # 7. DRY/WET DAY DIAGNOSTICS
  # ==========================================================================

  dryday_metrics <- compute_dryday_diagnostics(
    prcp_ref = prcp_ref,
    prcp_adj = prcp_adj,
    wet_thresh = wet_thresh
  )

  # ==========================================================================
  # 8. OVERALL SUMMARY
  # ==========================================================================

  summary_metrics <- compute_summary_metrics(
    moments = moments,
    quant_metrics = quant_metrics,
    extreme_metrics = extreme_metrics,
    spell_metrics = spell_metrics,
    dryday_metrics = dryday_metrics
  )

  # ==========================================================================
  # RETURN RESULTS
  # ==========================================================================

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
#' Computes validation diagnostics comparing original and adjusted precipitation
#' series using the same metrics returned by \code{\link{diagnose_prcp_qm}}.
#'
#' @details
#' This is a convenience wrapper around \code{\link{diagnose_prcp_qm}} with
#' argument names aligned to quantile-mapping workflows.
#'
#' @param prcp_org Numeric vector. Original precipitation values.
#' @param prcp_adjusted Numeric vector. Adjusted precipitation values (same length as \code{prcp_org}).
#' @param month Integer vector (optional). Month index (1--12) for each observation.
#' @param year Integer vector (optional). Year index (1..n_years) for each observation.
#' @param mean_factor Numeric matrix (optional). Target mean scaling factors (n_years x 12).
#' @param var_factor Numeric matrix (optional). Target variance scaling factors (n_years x 12).
#' @param wet_thresh Numeric. Threshold for defining wet days (default = 0.1 mm).
#' @param probs Numeric vector. Quantile probabilities to evaluate.
#'
#' @return A list of class \code{"prcp_qm_diagnostics"} with diagnostic tables
#'   and summary information (see \code{\link{diagnose_prcp_qm}} for details).
#'
#' @examples
#' prcp_org <- c(0, 1, 2, 0, 5, 3)
#' prcp_adjusted <- c(0, 1.1, 2.2, 0, 4.8, 3.1)
#' month <- c(1, 1, 1, 1, 1, 1)
#' year <- c(1, 1, 1, 1, 1, 1)
#' validate_quantile_mapping(prcp_org, prcp_adjusted, month = month, year = year)
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
# HELPER FUNCTIONS FOR DIAGNOSTICS
# ==============================================================================

#' Compute moment diagnostics for quantile-mapped precipitation
#' @keywords internal
compute_moment_diagnostics <- function(prcp_ref, prcp_adj, mask) {

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
  overall$diff <- overall$adjusted - overall$original
  overall$pct_change <- (overall$ratio - 1) * 100

  overall$assessment <- assess_moment_changes(overall)

  overall
}


#' Compute quantile diagnostics for quantile-mapped precipitation
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
  result$diff <- result$adjusted - result$original
  result$pct_change <- (result$ratio - 1) * 100
  result$abs_error <- abs(result$diff)

  result$assessment <- ifelse(
    abs(result$pct_change) < 5, "good",
    ifelse(abs(result$pct_change) < 15, "acceptable", "poor")
  )

  result
}


#' Compute extreme-value diagnostics for quantile-mapped precipitation
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


#' Compute temporal diagnostics for quantile-mapped precipitation
#' @keywords internal
compute_temporal_diagnostics <- function(prcp_ref, prcp_adj, month, year, mask) {

  annual_ref <- tapply(prcp_ref[mask], year[mask], mean, na.rm = TRUE)
  annual_adj <- tapply(prcp_adj[mask], year[mask], mean, na.rm = TRUE)

  monthly_ref <- tapply(prcp_ref[mask], month[mask], mean, na.rm = TRUE)
  monthly_adj <- tapply(prcp_adj[mask], month[mask], mean, na.rm = TRUE)

  cor_annual <- stats::cor(annual_ref, annual_adj, use = "complete.obs")
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
    assessment = c(
      ifelse(cor_annual > 0.9, "excellent", ifelse(cor_annual > 0.7, "good", "poor")),
      ifelse(cor_monthly > 0.95, "excellent", ifelse(cor_monthly > 0.85, "good", "poor")),
      NA, NA,
      ifelse(abs(acf_adj - acf_ref) < 0.1, "good", "moderate")
    ),
    stringsAsFactors = FALSE
  )
}


#' Compute month-by-month diagnostics for quantile-mapped precipitation
#' @keywords internal
compute_monthly_diagnostics <- function(prcp_ref, prcp_adj, month,
                                        mean_factor = NULL, var_factor = NULL,
                                        year = NULL) {

  monthly <- data.frame(
    month = 1:12,
    stringsAsFactors = FALSE
  )

  for (m in 1:12) {
    idx <- month == m & prcp_ref > 0 & prcp_adj > 0

    if (sum(idx) > 0) {
      monthly$original_mean[m] <- mean(prcp_ref[idx], na.rm = TRUE)
      monthly$adjusted_mean[m] <- mean(prcp_adj[idx], na.rm = TRUE)
      monthly$original_var[m] <- stats::var(prcp_ref[idx], na.rm = TRUE)
      monthly$adjusted_var[m] <- stats::var(prcp_adj[idx], na.rm = TRUE)
      monthly$n_obs[m] <- sum(idx)
    } else {
      monthly$original_mean[m] <- NA_real_
      monthly$adjusted_mean[m] <- NA_real_
      monthly$original_var[m] <- NA_real_
      monthly$adjusted_var[m] <- NA_real_
      monthly$n_obs[m] <- 0
    }
  }

  monthly$mean_ratio <- monthly$adjusted_mean / monthly$original_mean
  monthly$var_ratio <- monthly$adjusted_var / monthly$original_var

  if (!is.null(mean_factor) && !is.null(var_factor) && !is.null(year)) {
    monthly$target_mean_ratio <- colMeans(mean_factor, na.rm = TRUE)
    monthly$target_var_ratio <- colMeans(var_factor, na.rm = TRUE)

    monthly$mean_error <- abs(monthly$mean_ratio - monthly$target_mean_ratio)
    monthly$var_error <- abs(monthly$var_ratio - monthly$target_var_ratio)

    monthly$mean_assessment <- ifelse(
      monthly$mean_error < 0.05, "excellent",
      ifelse(monthly$mean_error < 0.15, "good", "poor")
    )

    monthly$var_assessment <- ifelse(
      monthly$var_error < 0.10, "excellent",
      ifelse(monthly$var_error < 0.25, "good", "poor")
    )
  }

  monthly
}


#' Compute wet and dry spell-length diagnostics
#' @keywords internal
compute_spell_diagnostics <- function(prcp_ref, prcp_adj, wet_thresh) {

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
  result$max_ratio <- result$adjusted_max / result$original_max
  result$sd_ratio <- result$adjusted_sd / result$original_sd

  result$assessment <- ifelse(
    abs(result$mean_ratio - 1) < 0.2, "good",
    ifelse(abs(result$mean_ratio - 1) < 0.4, "acceptable", "poor")
  )

  result
}


#' Compute dry-day and wet-day frequency diagnostics
#' @keywords internal
compute_dryday_diagnostics <- function(prcp_ref, prcp_adj, wet_thresh) {

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

  result$difference <- result$adjusted - result$original
  result$pct_change <- (result$adjusted / result$original - 1) * 100

  result$assessment <- c(
    ifelse(result$difference[1] == 0, "perfect", "changed"),
    ifelse(result$difference[2] == 0, "perfect", "changed"),
    ifelse(abs(result$difference[3]) < 0.01, "excellent", "changed"),
    ifelse(abs(result$difference[4]) < 0.01, "excellent", "changed")
  )

  result
}


#' Compute overall summary metrics for quantile-mapping diagnostics
#' @keywords internal
compute_summary_metrics <- function(moments, quant_metrics, extreme_metrics,
                                    spell_metrics, dryday_metrics) {

  out <- list()
  scores <- numeric()

  mean_score <- ifelse(abs(moments$pct_change[moments$metric == "mean"]) < 5, 10,
                       ifelse(abs(moments$pct_change[moments$metric == "mean"]) < 15, 5, 0))
  var_score <- ifelse(abs(moments$pct_change[moments$metric == "variance"]) < 10, 10,
                      ifelse(abs(moments$pct_change[moments$metric == "variance"]) < 25, 5, 0))
  cv_score <- ifelse(abs(moments$pct_change[moments$metric == "cv"]) < 10, 10,
                     ifelse(abs(moments$pct_change[moments$metric == "cv"]) < 20, 5, 0))
  scores["moments"] <- mean_score + var_score + cv_score

  good_q <- sum(quant_metrics$assessment == "good")
  ok_q <- sum(quant_metrics$assessment == "acceptable")
  q_score <- (good_q * 4 + ok_q * 2)
  scores["quantiles"] <- min(q_score, 30)

  ext_score <- 20 * (1 - mean(abs(extreme_metrics$mean_ratio - 1)))
  scores["extremes"] <- max(0, min(ext_score, 20))

  sp_score <- sum(spell_metrics$assessment == "good") * 5
  scores["spells"] <- min(sp_score, 10)

  dry_score <- ifelse(dryday_metrics$assessment[dryday_metrics$category == "dry_days"] == "perfect", 10, 0)
  scores["drydays"] <- dry_score

  out$overall_score <- sum(scores)
  out$component_scores <- scores

  out$overall_assessment <- ifelse(
    out$overall_score >= 80, "excellent",
    ifelse(out$overall_score >= 60, "good",
           ifelse(out$overall_score >= 40, "acceptable", "poor"))
  )

  out$key_findings <- list(
    mean_change_pct = moments$pct_change[moments$metric == "mean"],
    var_change_pct = moments$pct_change[moments$metric == "variance"],
    dryday_preserved = dryday_metrics$difference[dryday_metrics$category == "dry_days"] == 0,
    worst_prob = quant_metrics$prob[which.max(abs(quant_metrics$pct_change))],
    extreme_amplification = mean(extreme_metrics$mean_ratio)
  )

  out
}


# ==============================================================================
# PRINT AND PLOT METHODS
# ==============================================================================

#' Print precipitation QM diagnostics
#'
#' @description
#' Prints a compact summary of quantile-mapping diagnostics.
#'
#' @param x A \code{"prcp_qm_diagnostics"} object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.prcp_qm_diagnostics <- function(x, ...) {
  cat("\n=== Precipitation QM Diagnostics ===\n\n")

  cat("Overall Assessment:", x$summary$overall_assessment, "\n")
  cat("Overall Score:", round(x$summary$overall_score, 1), "/ 100\n\n")

  cat("Component Scores:\n")
  for (comp in names(x$summary$component_scores)) {
    cat(sprintf("  %s: %.1f\n", comp, x$summary$component_scores[comp]))
  }
  cat("\n")

  cat("Key Findings:\n")
  cat(sprintf("  Mean change: %.1f%%\n", x$summary$key_findings$mean_change_pct))
  cat(sprintf("  Variance change: %.1f%%\n", x$summary$key_findings$var_change_pct))
  cat(sprintf("  Dry days preserved: %s\n", x$summary$key_findings$dryday_preserved))
  cat(sprintf("  Extreme amplification: %.2fx\n", x$summary$key_findings$extreme_amplification))
  cat("\n")

  cat("Moment Preservation:\n")
  print(x$moments[, c("metric", "original", "adjusted", "pct_change", "assessment")],
        row.names = FALSE, digits = 3)
  cat("\n")

  cat("Dry/Wet Day Preservation:\n")
  print(x$drydays[x$drydays$category %in% c("dry_days", "wet_days"),
                  c("category", "original", "adjusted", "assessment")],
        row.names = FALSE)
  cat("\n")

  invisible(x)
}


#' Plot diagnostics for precipitation QM
#'
#' @description
#' Plots selected diagnostic panels from a \code{"prcp_qm_diagnostics"} object.
#'
#' @param x A \code{"prcp_qm_diagnostics"} object.
#' @param which Character. Which plot to return.
#' @param ... Additional arguments passed to plotting functions (unused).
#'
#' @return A ggplot object when a single plot is requested, or a list of plots
#'   when \code{which = "all"} and \code{gridExtra} is unavailable.
#'
#' @export
plot.prcp_qm_diagnostics <- function(x, which = c("all", "moments", "quantiles",
                                                  "extremes", "monthly"), ...) {

  which <- match.arg(which)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }

  plots <- list()

  if (which %in% c("all", "moments")) {
    p1 <- ggplot(x$moments, aes(x = rlang::.data$metric)) +
      geom_col(aes(y = rlang::.data$pct_change, fill = rlang::.data$assessment)) +
      scale_fill_manual(
        values = c("excellent" = "green3", "good" = "yellow3", "poor" = "red3")
      ) +
      labs(
        title = "Moment Changes",
        x = "Metric",
        y = "Percent Change (%)",
        fill = "Assessment"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plots$moments <- p1
  }

  if (which %in% c("all", "quantiles")) {
    p2 <- ggplot(x$quantiles, aes(x = original, y = adjusted)) +
      geom_point(aes(color = assessment), size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      scale_color_manual(
        values = c("good" = "green3", "acceptable" = "yellow3", "poor" = "red3")
      ) +
      labs(
        title = "Quantile Preservation",
        x = "Reference quantile value",
        y = "Adjusted quantile value",
        color = "Assessment"
      ) +
      theme_minimal()

    plots$quantiles <- p2
  }

  if (which %in% c("all", "extremes")) {
    extreme_long <- data.frame(
      threshold = rep(x$extremes$threshold, 2),
      type = rep(c("Reference", "Adjusted"), each = nrow(x$extremes)),
      mean = c(x$extremes$original_mean, x$extremes$adjusted_mean)
    )

    p3 <- ggplot(extreme_long, aes(
      x = factor(threshold),
      y = mean,
      fill = type
    )) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("Reference" = "blue3", "Adjusted" = "red3")) +
      labs(
        title = "Extreme Value Comparison",
        x = "Exceedance quantile",
        y = "Mean value above cutoff",
        fill = ""
      ) +
      theme_minimal()

    plots$extremes <- p3
  }

  if (which %in% c("all", "monthly") && !is.null(x$monthly)) {
    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("Package 'tidyr' required for monthly plot")
    }

    monthly_long <- tidyr::pivot_longer(
      x$monthly[, c("month", "original_mean", "adjusted_mean")],
      cols = c("original_mean", "adjusted_mean"),
      names_to = "type",
      values_to = "mean"
    )
    monthly_long$type <- ifelse(monthly_long$type == "original_mean", "Reference", "Adjusted")

    p4 <- ggplot(monthly_long, aes(
      x = factor(month),
      y = mean,
      color = type,
      group = type
    )) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_color_manual(values = c("Reference" = "blue3", "Adjusted" = "red3")) +
      scale_x_discrete(labels = month.abb) +
      labs(
        title = "Monthly Mean Precipitation",
        x = "Month",
        y = "Mean (mm/day)",
        color = ""
      ) +
      theme_minimal()

    plots$monthly <- p4
  }

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


#' Summarize precipitation QM diagnostics
#'
#' @description
#' Prints a text summary for a \code{"prcp_qm_diagnostics"} object.
#'
#' @param object A \code{"prcp_qm_diagnostics"} object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
summary.prcp_qm_diagnostics <- function(object, ...) {
  cat("\n=== Precipitation QM Diagnostics Summary ===\n\n")

  cat("Overall Score:", round(object$summary$overall_score, 1), "/ 100\n")
  cat("Assessment:", object$summary$overall_assessment, "\n\n")

  cat("Detailed Assessment:\n")
  cat("------------------\n")

  cat("\nMoments:\n")
  good_moments <- sum(object$moments$assessment %in% c("excellent", "good"))
  cat(sprintf(
    "  %s/%s metrics within acceptable range\n",
    format(good_moments, big.mark = ","),
    format(nrow(object$moments), big.mark = ",")
  ))

  cat("\nQuantiles:\n")
  good_quant <- sum(object$quantiles$assessment == "good")
  cat(sprintf(
    "  %s/%s quantiles well preserved\n",
    format(good_quant, big.mark = ","),
    format(nrow(object$quantiles), big.mark = ",")
  ))

  cat("\nDry Days:\n")
  dd <- object$drydays
  dd_assess <- dd$assessment[dd$category == "dry_days"]
  cat(sprintf("  Preservation: %s\n", dd_assess))

  cat("\nExtremes:\n")
  cat(sprintf("  Mean amplification: %.2fx\n", object$summary$key_findings$extreme_amplification))

  invisible(object)
}
