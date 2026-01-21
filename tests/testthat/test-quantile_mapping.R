# ---- adjust_precipitation_qm -------------------------------------------------

testthat::test_that("adjust_precipitation_qm: year index contiguity is enforced (guards against calendar years)", {
  set.seed(10)

  n_years <- 2L
  n_days <- 365L * n_years
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.2, scale = 5)
  precip[sample.int(n_days, size = 80)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  # Non-contiguous year indices (mimics passing calendar years)
  year_bad <- rep(c(2020L, 2021L), each = 365L)

  testthat::expect_error(
    adjust_precipitation_qm(
      precip = precip,
      mean_factor = mean_factor,
      var_factor = var_factor,
      month = month,
      year = year_bad,
      verbose = FALSE
    ),
    "contiguous simulation-year index"
  )
})

testthat::test_that("adjust_precipitation_qm: identity mapping (no mean enforcement) preserves wet days and dry days", {
  set.seed(1)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.5, scale = 4)
  precip[sample.int(n_days, size = 120)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    intensity_threshold = 0,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    min_events = 10,
    validate_output = TRUE,
    diagnostics = FALSE,
    verbose = FALSE
  )

  testthat::expect_type(out, "double")
  testthat::expect_length(out, length(precip))
  testthat::expect_true(all(is.finite(out[!is.na(out)])))

  # Dry days remain unchanged
  testthat::expect_identical(out[precip == 0], precip[precip == 0])

  # With identical params and no enforcement, qgamma(pgamma(x)) should be ~ x
  wet <- precip > 0
  testthat::expect_lt(mean(abs(out[wet] - precip[wet])), 1e-6)

  # Attributes exist
  testthat::expect_true(!is.null(attr(out, "perturbed_months")))
  testthat::expect_true(!is.null(attr(out, "skipped_months")))
  testthat::expect_true(!is.null(attr(out, "n_failed_fits")))
})

testthat::test_that("adjust_precipitation_qm: mean_factor scales wet-day mean approximately (no tail exaggeration)", {
  set.seed(2)

  n_years <- 3L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.3, scale = 5)
  precip[sample.int(n_days, size = 150)] <- 0

  mean_factor <- matrix(0.7, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    intensity_threshold = 0,
    min_events = 10,
    validate_output = TRUE,
    verbose = FALSE
  )

  # Check mean scaling on wet days across all months/years combined
  wet <- precip > 0
  m0 <- mean(precip[wet])
  m1 <- mean(out[wet])

  # Allow tolerance because mapping is month-wise and depends on fits
  testthat::expect_equal(m1 / m0, 0.7, tolerance = 0.05)

  # Dry days unchanged
  testthat::expect_identical(out[precip == 0], precip[precip == 0])
})

testthat::test_that("adjust_precipitation_qm: variance_factor affects spread of wet-day intensities", {
  set.seed(3)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.2, scale = 6)
  precip[sample.int(n_days, size = 100)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.6, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    min_events = 10,
    validate_output = TRUE,
    verbose = FALSE
  )

  wet <- precip > 0
  v0 <- stats::var(precip[wet])
  v1 <- stats::var(out[wet])

  # Should increase variance noticeably (not necessarily exactly by factor, due to month-wise mapping)
  testthat::expect_gt(v1, v0 * 1.2)
})

testthat::test_that("adjust_precipitation_qm: scale_var_with_mean combines variance and mean scaling", {
  set.seed(4)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.8, scale = 3)
  precip[sample.int(n_days, size = 80)] <- 0

  mean_factor <- matrix(0.8, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out_noscale <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    scale_var_with_mean = FALSE,
    month = month,
    year = year,
    verbose = FALSE
  )

  out_scale <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    scale_var_with_mean = TRUE,
    month = month,
    year = year,
    verbose = FALSE
  )

  wet <- precip > 0

  # Means should be similar (both target mean ~ 0.8 * baseline)
  testthat::expect_equal(mean(out_scale[wet]) / mean(precip[wet]), 0.8, tolerance = 0.05)
  testthat::expect_equal(mean(out_noscale[wet]) / mean(precip[wet]), 0.8, tolerance = 0.05)

  # With scale_var_with_mean=TRUE, effective variance factor ~ mean_factor^2,
  # so variance should be lower than the no-scale case (which keeps variance factor at 1).
  testthat::expect_lt(stats::var(out_scale[wet]), stats::var(out_noscale[wet]))
})

testthat::test_that("adjust_precipitation_qm: tail amplification increases upper tail; mean enforcement keeps mean on target", {
  set.seed(5)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.1, scale = 7)
  precip[sample.int(n_days, size = 90)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out_base <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    exaggerate_extremes = FALSE,
    month = month,
    year = year,
    verbose = FALSE
  )

  out_tail_nom <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    exaggerate_extremes = TRUE,
    extreme_prob_threshold = 0.95,
    extreme_k = 1.4,
    enforce_target_mean = FALSE,
    month = month,
    year = year,
    verbose = FALSE
  )

  out_tail_mean <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    exaggerate_extremes = TRUE,
    extreme_prob_threshold = 0.95,
    extreme_k = 1.4,
    enforce_target_mean = TRUE,
    month = month,
    year = year,
    verbose = FALSE
  )

  wet <- precip > 0

  q0 <- stats::quantile(out_base[wet], probs = 0.99, names = FALSE, na.rm = TRUE)
  q1 <- stats::quantile(out_tail_nom[wet], probs = 0.99, names = FALSE, na.rm = TRUE)
  q2 <- stats::quantile(out_tail_mean[wet], probs = 0.99, names = FALSE, na.rm = TRUE)

  # Tail amplification should increase the upper tail relative to baseline
  testthat::expect_gt(q1, q0)
  testthat::expect_gt(q2, q0)

  # If enforce_target_mean=TRUE and mean_factor=1, wet-day mean should remain ~ baseline wet-day mean
  testthat::expect_equal(mean(out_tail_mean[wet]), mean(out_base[wet]), tolerance = 0.02)

  # Without mean enforcement, mean can drift (not guaranteed, but expected with tail exaggeration)
  # Require that drift is non-zero beyond tiny numerical error.
  testthat::expect_true(abs(mean(out_tail_nom[wet]) - mean(out_base[wet])) > 1e-6)
})

testthat::test_that("adjust_precipitation_qm: intensity_threshold keeps small values unchanged and does not change dry-day frequency", {
  set.seed(6)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.4, scale = 4)

  # Inject a band of small drizzle values that should be treated as "dry" for intensity mapping
  drizzle_idx <- sample.int(n_days, size = 120)
  precip[drizzle_idx] <- runif(length(drizzle_idx), min = 0, max = 0.5)

  # Inject true dry days too
  precip[sample(setdiff(seq_len(n_days), drizzle_idx), size = 60)] <- 0

  thr <- 0.5

  mean_factor <- matrix(1.2, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    intensity_threshold = thr,
    verbose = FALSE
  )

  # Values <= threshold must be unchanged (including exact zeros and drizzle band)
  keep <- !is.na(precip) & (precip <= thr)
  testthat::expect_identical(out[keep], precip[keep])

  # Wet-day indicator relative to threshold is preserved (function does not change occurrence)
  wet0 <- precip > thr
  wet1 <- out > thr
  wet0[is.na(wet0)] <- FALSE
  wet1[is.na(wet1)] <- FALSE
  testthat::expect_identical(wet1, wet0)
})

testthat::test_that("adjust_precipitation_qm: months with insufficient events are skipped and pass through", {
  set.seed(7)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)

  # Force month=1 to have very few wet days by setting nearly all to 0 in month 1
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.6, scale = 4)
  precip[sample.int(n_days, size = 120)] <- 0

  # Make month 1 mostly dry, so it will fail min_events for wet intensities
  idx_m1 <- which(month == 1 & precip > 0)
  if (length(idx_m1) > 5) precip[idx_m1[-seq_len(5)]] <- 0

  mean_factor <- matrix(1.3, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    min_events = 10,
    verbose = FALSE
  )

  # Any wet days in skipped month should pass through unchanged
  wet_m1 <- which(month == 1 & precip > 0)
  if (length(wet_m1) > 0) {
    testthat::expect_identical(out[wet_m1], precip[wet_m1])
  }

  skipped <- attr(out, "skipped_months")
  testthat::expect_true(1 %in% skipped)
})

testthat::test_that("adjust_precipitation_qm: NA values pass through unchanged", {
  set.seed(8)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.3, scale = 5)
  precip[sample.int(n_days, size = 80)] <- 0
  na_idx <- sample.int(n_days, size = 20)
  precip[na_idx] <- NA_real_

  mean_factor <- matrix(1.1, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    verbose = FALSE
  )

  testthat::expect_true(all(is.na(out[na_idx])))
})

testthat::test_that("adjust_precipitation_qm: diagnostics=TRUE returns expected structure", {
  set.seed(9)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.4, scale = 4)
  precip[sample.int(n_days, size = 100)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  res <- adjust_precipitation_qm(
    precip = precip,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    diagnostics = TRUE,
    verbose = FALSE
  )

  testthat::expect_true(is.list(res))
  testthat::expect_true(all(c("adjusted", "diagnostics") %in% names(res)))

  out <- res$adjusted
  testthat::expect_length(out, length(precip))
  testthat::expect_true(inherits(res$diagnostics, "precip_qm_diagnostics"))

  # Attributes exist on adjusted output
  testthat::expect_true(!is.null(attr(out, "perturbed_months")))
  testthat::expect_true(!is.null(attr(out, "skipped_months")))
  testthat::expect_true(!is.null(attr(out, "n_failed_fits")))
})

# ---- diagnose_precip_qm / validate_quantile_mapping --------------------------

testthat::test_that("diagnose_precip_qm: returns expected structure and intended ratios", {
  set.seed(11)

  n_years <- 2L
  n_days <- 12L * n_years
  month <- rep(1:12, n_years)
  year <- rep(seq_len(n_years), each = 12L)

  precip_ref <- rgamma(n_days, shape = 2, scale = 2) + 0.2
  precip_ref[c(3, 14)] <- 0
  precip_adj <- precip_ref * 1.2
  precip_adj[precip_ref == 0] <- 0

  mean_factor <- matrix(1.2, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.44, nrow = n_years, ncol = 12)

  diag <- diagnose_precip_qm(
    precip_ref = precip_ref,
    precip_adj = precip_adj,
    month = month,
    year = year,
    mean_factor = mean_factor,
    var_factor = var_factor,
    wet_thresh = 0.1,
    probs = c(0.5, 0.9)
  )

  testthat::expect_s3_class(diag, "precip_qm_diagnostics")
  testthat::expect_true(all(c("moments", "quantiles", "extremes", "spells", "drydays", "summary") %in% names(diag)))
  testthat::expect_s3_class(diag$moments, "data.frame")
  testthat::expect_s3_class(diag$monthly, "data.frame")
  testthat::expect_true(nrow(diag$monthly) == 12L)

  mean_row <- diag$moments[diag$moments$metric == "mean", , drop = FALSE]
  var_row  <- diag$moments[diag$moments$metric == "variance", , drop = FALSE]

  testthat::expect_equal(mean_row$intended_ratio, 1.2, tolerance = 1e-8)
  testthat::expect_equal(var_row$intended_ratio, 1.44, tolerance = 1e-8)
  testthat::expect_equal(mean_row$ratio, 1.2, tolerance = 1e-6)
  testthat::expect_equal(var_row$ratio, 1.44, tolerance = 1e-6)

  dry_days <- diag$drydays[diag$drydays$category == "dry_days", , drop = FALSE]
  testthat::expect_identical(dry_days$diff, 0)
})

testthat::test_that("validate_quantile_mapping: matches diagnose_precip_qm output and checks length", {
  set.seed(12)

  n_years <- 2L
  n_days <- 12L * n_years
  month <- rep(1:12, n_years)
  year <- rep(seq_len(n_years), each = 12L)

  precip_ref <- rgamma(n_days, shape = 2.1, scale = 2) + 0.2
  precip_ref[c(2, 18)] <- 0
  precip_adj <- precip_ref * 0.9
  precip_adj[precip_ref == 0] <- 0

  mean_factor <- matrix(0.9, nrow = n_years, ncol = 12)
  var_factor  <- matrix(0.81, nrow = n_years, ncol = 12)

  diag_a <- diagnose_precip_qm(
    precip_ref = precip_ref,
    precip_adj = precip_adj,
    month = month,
    year = year,
    mean_factor = mean_factor,
    var_factor = var_factor,
    wet_thresh = 0.1,
    probs = c(0.5)
  )

  diag_b <- validate_quantile_mapping(
    precip_org = precip_ref,
    precip_adjusted = precip_adj,
    month = month,
    year = year,
    mean_factor = mean_factor,
    var_factor = var_factor,
    wet_thresh = 0.1,
    probs = c(0.5)
  )

  testthat::expect_s3_class(diag_b, "precip_qm_diagnostics")
  testthat::expect_equal(diag_b$summary, diag_a$summary)
  testthat::expect_equal(diag_b$quantiles$prob, diag_a$quantiles$prob)

  testthat::expect_error(
    validate_quantile_mapping(
      precip_org = 1:3,
      precip_adjusted = 1:2
    ),
    "must have same length"
  )
})
