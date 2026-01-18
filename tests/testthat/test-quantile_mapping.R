testthat::test_that("adjust_precipitation_qm: identity mapping (no mean enforcement) preserves wet days and dry days", {
  set.seed(1)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  prcp <- rgamma(n_days, shape = 1.5, scale = 4)
  prcp[sample.int(n_days, size = 120)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    intensity_threshold = 0,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,   # critical for identity behavior
    min_events = 10,
    validate_output = TRUE,
    diagnostics = FALSE,
    verbose = FALSE
  )

  testthat::expect_type(out, "double")
  testthat::expect_length(out, length(prcp))
  testthat::expect_true(all(is.finite(out[!is.na(out)])))

  # Dry days remain unchanged
  testthat::expect_identical(out[prcp == 0], prcp[prcp == 0])

  # With identical params and no enforcement, qgamma(pgamma(x)) should be ~ x
  wet <- prcp > 0
  testthat::expect_lt(mean(abs(out[wet] - prcp[wet])), 1e-6)

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

  prcp <- rgamma(n_days, shape = 1.3, scale = 5)
  prcp[sample.int(n_days, size = 150)] <- 0

  mean_factor <- matrix(0.7, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    prcp = prcp,
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
  wet <- prcp > 0
  m0 <- mean(prcp[wet])
  m1 <- mean(out[wet])

  # Allow tolerance because mapping is month-wise and depends on fits
  testthat::expect_equal(m1 / m0, 0.7, tolerance = 0.05)

  # Dry days unchanged
  testthat::expect_identical(out[prcp == 0], prcp[prcp == 0])
})

testthat::test_that("adjust_precipitation_qm: variance_factor affects spread of wet-day intensities", {
  set.seed(3)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  prcp <- rgamma(n_days, shape = 1.2, scale = 6)
  prcp[sample.int(n_days, size = 100)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.6, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    min_events = 10,
    validate_output = TRUE,
    verbose = FALSE
  )

  wet <- prcp > 0
  v0 <- stats::var(prcp[wet])
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

  prcp <- rgamma(n_days, shape = 1.8, scale = 3)
  prcp[sample.int(n_days, size = 80)] <- 0

  mean_factor <- matrix(0.8, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out_noscale <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    scale_var_with_mean = FALSE,
    month = month,
    year = year,
    verbose = FALSE
  )

  out_scale <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    scale_var_with_mean = TRUE,
    month = month,
    year = year,
    verbose = FALSE
  )

  wet <- prcp > 0

  # Means should be similar (both target mean ~ 0.8 * baseline)
  testthat::expect_equal(mean(out_scale[wet]) / mean(prcp[wet]), 0.8, tolerance = 0.05)
  testthat::expect_equal(mean(out_noscale[wet]) / mean(prcp[wet]), 0.8, tolerance = 0.05)

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

  prcp <- rgamma(n_days, shape = 1.1, scale = 7)
  prcp[sample.int(n_days, size = 90)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out_base <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    exaggerate_extremes = FALSE,
    month = month,
    year = year,
    verbose = FALSE
  )

  out_tail_nom <- adjust_precipitation_qm(
    prcp = prcp,
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
    prcp = prcp,
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

  wet <- prcp > 0

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

  prcp <- rgamma(n_days, shape = 1.4, scale = 4)

  # Inject a band of small drizzle values that should be treated as "dry" for intensity mapping
  drizzle_idx <- sample.int(n_days, size = 120)
  prcp[drizzle_idx] <- runif(length(drizzle_idx), min = 0, max = 0.5)

  # Inject true dry days too
  prcp[sample(setdiff(seq_len(n_days), drizzle_idx), size = 60)] <- 0

  thr <- 0.5

  mean_factor <- matrix(1.2, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    intensity_threshold = thr,
    verbose = FALSE
  )

  # Values <= threshold must be unchanged (including exact zeros and drizzle band)
  keep <- !is.na(prcp) & (prcp <= thr)
  testthat::expect_identical(out[keep], prcp[keep])

  # Wet-day indicator relative to threshold is preserved (function does not change occurrence)
  wet0 <- prcp > thr
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

  prcp <- rgamma(n_days, shape = 1.6, scale = 4)
  prcp[sample.int(n_days, size = 120)] <- 0

  # Make month 1 mostly dry, so it will fail min_events for wet intensities
  idx_m1 <- which(month == 1 & prcp > 0)
  if (length(idx_m1) > 5) prcp[idx_m1[-seq_len(5)]] <- 0

  mean_factor <- matrix(1.3, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    min_events = 10,
    verbose = FALSE
  )

  # Any wet days in skipped month should pass through unchanged
  wet_m1 <- which(month == 1 & prcp > 0)
  if (length(wet_m1) > 0) {
    testthat::expect_identical(out[wet_m1], prcp[wet_m1])
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

  prcp <- rgamma(n_days, shape = 1.3, scale = 5)
  prcp[sample.int(n_days, size = 80)] <- 0
  na_idx <- sample.int(n_days, size = 20)
  prcp[na_idx] <- NA_real_

  mean_factor <- matrix(1.1, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  out <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    verbose = FALSE
  )

  testthat::expect_true(all(is.na(out[na_idx])))
})

testthat::test_that("adjust_precipitation_qm: diagnostics=TRUE returns expected structure and includes fitted objects", {
  set.seed(9)

  n_years <- 2L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  prcp <- rgamma(n_days, shape = 1.4, scale = 4)
  prcp[sample.int(n_days, size = 100)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  res <- adjust_precipitation_qm(
    prcp = prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = month,
    year = year,
    diagnostics = TRUE,
    verbose = FALSE
  )

  testthat::expect_true(is.list(res))
  testthat::expect_true(all(c("adjusted", "diagnostics", "base_gamma", "target_gamma", "var_factor_use") %in% names(res)))

  out <- res$adjusted
  testthat::expect_length(out, length(prcp))

  # base_gamma is data.frame with required columns
  testthat::expect_s3_class(res$base_gamma, "data.frame")
  testthat::expect_true(all(c("month", "shape", "scale", "mean", "var") %in% names(res$base_gamma)))

  # target_gamma list contains matrices indexed [month_row, year_index]
  testthat::expect_true(is.list(res$target_gamma))
  testthat::expect_true(all(c("months", "shape", "scale", "mean", "var") %in% names(res$target_gamma)))
  testthat::expect_true(is.matrix(res$target_gamma$shape))
  testthat::expect_true(ncol(res$target_gamma$shape) == n_years)

  # var_factor_use same dims as inputs
  testthat::expect_true(is.matrix(res$var_factor_use))
  testthat::expect_true(all(dim(res$var_factor_use) == c(n_years, 12)))
})

testthat::test_that("adjust_precipitation_qm: year index contiguity is enforced (guards against calendar years)", {
  set.seed(10)

  n_years <- 2L
  n_days <- 365L * n_years
  month <- rep(1:12, length.out = n_days)

  prcp <- rgamma(n_days, shape = 1.2, scale = 5)
  prcp[sample.int(n_days, size = 80)] <- 0

  mean_factor <- matrix(1.0, nrow = n_years, ncol = 12)
  var_factor  <- matrix(1.0, nrow = n_years, ncol = 12)

  # Non-contiguous year indices (mimics passing calendar years)
  year_bad <- rep(c(2020L, 2021L), each = 365L)

  testthat::expect_error(
    adjust_precipitation_qm(
      prcp = prcp,
      mean_factor = mean_factor,
      var_factor = var_factor,
      month = month,
      year = year_bad,
      verbose = FALSE
    ),
    "contiguous simulation-year index"
  )
})
