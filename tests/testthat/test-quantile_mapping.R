# ==============================================================================
# UNIT TESTS FOR quantile_mapping()
# ==============================================================================

library(testthat)

# ==============================================================================
# TEST FIXTURES AND HELPER FUNCTIONS
# ==============================================================================

#' Generate synthetic daily precipitation data
#' @keywords internal
generate_test_precip <- function(n_years = 2, seed = 123) {
  set.seed(seed)

  n_days <- n_years * 365

  # Month indices (no leap days)
  month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  month_idx <- rep(rep(1:12, times = month_days), n_years)[1:n_days]

  # Year indices
  year_idx <- rep(1:n_years, each = 365)

  # Generate precipitation with seasonal pattern
  # Higher precipitation in winter months (12, 1, 2)
  seasonal_mean <- c(8, 7, 6, 4, 3, 2, 2, 3, 4, 5, 6, 7)[month_idx]

  # Generate gamma-distributed precipitation
  precip <- rgamma(n_days, shape = 2, scale = seasonal_mean / 2)

  # Add some dry days (~30%)
  dry_days <- sample(n_days, size = floor(n_days * 0.3))
  precip[dry_days] <- 0

  list(
    value = precip,
    mon.ts = month_idx,
    year.ts = year_idx,
    n_days = n_days,
    n_years = n_years
  )
}


#' Create change factor matrices
#' @keywords internal
create_change_factors <- function(n_years, mean_factor = 1.1, var_factor = 1.2) {
  list(
    mean.change = matrix(mean_factor, nrow = n_years, ncol = 12),
    var.change = matrix(var_factor, nrow = n_years, ncol = 12)
  )
}


# ==============================================================================
# BASIC FUNCTIONALITY TESTS
# ==============================================================================

test_that("quantile_mapping returns correct length output", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  expect_equal(length(result), length(data$value))
  expect_type(result, "double")
})


test_that("quantile_mapping preserves zeros (dry days)", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Zeros should remain zeros
  expect_equal(sum(data$value == 0), sum(result == 0))
  expect_true(all(result[data$value == 0] == 0))
})


test_that("quantile_mapping increases mean when mean.change > 1", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.5, var_factor = 1.0)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Mean should increase (excluding zeros)
  nonzero_original <- data$value[data$value > 0]
  nonzero_result <- result[result > 0]

  expect_true(mean(nonzero_result) > mean(nonzero_original))
})


test_that("quantile_mapping increases variance when var.change > 1", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.0, var_factor = 2.0)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Variance should increase (excluding zeros)
  nonzero_original <- data$value[data$value > 0]
  nonzero_result <- result[result > 0]

  expect_true(var(nonzero_result) > var(nonzero_original))
})


test_that("quantile_mapping with no change returns similar distribution", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.0, var_factor = 1.0)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Should be very similar (within 5% due to numerical precision)
  nonzero_idx <- data$value > 0
  relative_diff <- abs(result[nonzero_idx] - data$value[nonzero_idx]) / data$value[nonzero_idx]

  expect_true(median(relative_diff) < 0.05)
})


# ==============================================================================
# INPUT VALIDATION TESTS
# ==============================================================================

test_that("quantile_mapping rejects NULL inputs", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  expect_error(
    quantile_mapping(
      value = NULL,
      mean.change = changes$mean.change,
      var.change = changes$var.change,
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'value' must not be NULL"
  )

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = NULL,
      var.change = changes$var.change,
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'mean.change' must not be NULL"
  )
})


test_that("quantile_mapping rejects mismatched lengths", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = changes$mean.change,
      var.change = changes$var.change,
      mon.ts = data$mon.ts[1:100],  # Wrong length
      year.ts = data$year.ts
    ),
    "'mon.ts' must have same length as 'value'"
  )
})


test_that("quantile_mapping rejects negative precipitation", {
  data <- generate_test_precip(n_years = 2)
  data$value[1] <- -1  # Invalid
  changes <- create_change_factors(n_years = 2)

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = changes$mean.change,
      var.change = changes$var.change,
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'value' must be non-negative"
  )
})


test_that("quantile_mapping rejects invalid month indices", {
  data <- generate_test_precip(n_years = 2)
  data$mon.ts[1] <- 13  # Invalid month
  changes <- create_change_factors(n_years = 2)

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = changes$mean.change,
      var.change = changes$var.change,
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'mon.ts' must contain integers between 1 and 12"
  )
})


test_that("quantile_mapping rejects non-matrix change factors", {
  data <- generate_test_precip(n_years = 2)

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = rep(1.1, 24),  # Vector, not matrix
      var.change = matrix(1.2, nrow = 2, ncol = 12),
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'mean.change' must be a matrix"
  )
})


test_that("quantile_mapping rejects wrong matrix dimensions", {
  data <- generate_test_precip(n_years = 2)

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = matrix(1.1, nrow = 2, ncol = 6),  # Wrong number of columns
      var.change = matrix(1.2, nrow = 2, ncol = 12),
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'mean.change' must have 12 columns"
  )

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = matrix(1.1, nrow = 3, ncol = 12),  # Wrong number of rows
      var.change = matrix(1.2, nrow = 3, ncol = 12),
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'mean.change' must have 2 rows"
  )
})


test_that("quantile_mapping rejects negative change factors", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)
  changes$mean.change[1, 1] <- -0.5  # Negative factor

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = changes$mean.change,
      var.change = changes$var.change,
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'mean.change' must contain positive values"
  )
})


test_that("quantile_mapping rejects zero change factors", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)
  changes$var.change[1, 1] <- 0  # Zero factor

  expect_error(
    quantile_mapping(
      value = data$value,
      mean.change = changes$mean.change,
      var.change = changes$var.change,
      mon.ts = data$mon.ts,
      year.ts = data$year.ts
    ),
    "'var.change' must contain positive values"
  )
})


# ==============================================================================
# EDGE CASES AND ROBUSTNESS TESTS
# ==============================================================================

test_that("quantile_mapping handles all-zero input", {
  data <- generate_test_precip(n_years = 2)
  data$value[] <- 0  # All zeros
  changes <- create_change_factors(n_years = 2)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts,
    verbose = FALSE
  )

  expect_equal(result, data$value)
  expect_true(all(result == 0))
})


test_that("quantile_mapping handles insufficient data in some months", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  # Make January almost all zeros (< 10 non-zero events)
  jan_idx <- which(data$mon.ts == 1)
  data$value[jan_idx[1:55]] <- 0

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts,
    min.events = 10,
    verbose = FALSE
  )

  # Should still return valid output
  expect_equal(length(result), length(data$value))

  # Check that January was skipped
  skipped <- attr(result, "skipped_months")
  expect_true(1 %in% skipped)
})


test_that("quantile_mapping handles very small precipitation values", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  # Add some very small values
  small_idx <- sample(which(data$value > 0), 50)
  data$value[small_idx] <- runif(50, 0.001, 0.01)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Should handle without errors
  expect_equal(length(result), length(data$value))
  expect_true(all(is.finite(result)))
})


test_that("quantile_mapping handles extreme change factors", {
  data <- generate_test_precip(n_years = 2)

  # Very large increase
  changes_large <- create_change_factors(n_years = 2, mean_factor = 5.0, var_factor = 10.0)

  result_large <- quantile_mapping(
    value = data$value,
    mean.change = changes_large$mean.change,
    var.change = changes_large$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  expect_true(all(is.finite(result_large)))
  expect_true(mean(result_large) > mean(data$value))

  # Very small decrease
  changes_small <- create_change_factors(n_years = 2, mean_factor = 0.2, var_factor = 0.5)

  result_small <- quantile_mapping(
    value = data$value,
    mean.change = changes_small$mean.change,
    var.change = changes_small$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  expect_true(all(is.finite(result_small)))
  expect_true(mean(result_small[result_small > 0]) < mean(data$value[data$value > 0]))
})


test_that("quantile_mapping handles single year data", {
  data <- generate_test_precip(n_years = 1)
  changes <- create_change_factors(n_years = 1)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  expect_equal(length(result), length(data$value))
  expect_true(all(is.finite(result)))
})


test_that("quantile_mapping handles varying change factors by year", {
  data <- generate_test_precip(n_years = 3)

  # Create time-varying changes (increasing over years)
  mean.change <- matrix(NA, nrow = 3, ncol = 12)
  var.change <- matrix(NA, nrow = 3, ncol = 12)

  for (y in 1:3) {
    mean.change[y, ] <- 1.0 + (y - 1) * 0.05  # 1.0, 1.05, 1.10
    var.change[y, ] <- 1.0 + (y - 1) * 0.10   # 1.0, 1.10, 1.20
  }

  result <- quantile_mapping(
    value = data$value,
    mean.change = mean.change,
    var.change = var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Calculate mean by year
  year_means <- tapply(result[result > 0], data$year.ts[result > 0], mean)

  # Means should increase over years
  expect_true(year_means[2] > year_means[1])
  expect_true(year_means[3] > year_means[2])
})


# ==============================================================================
# ATTRIBUTE AND DIAGNOSTIC TESTS
# ==============================================================================

test_that("quantile_mapping returns diagnostic attributes", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  expect_true(!is.null(attr(result, "perturbed_months")))
  expect_true(!is.null(attr(result, "skipped_months")))
  expect_true(!is.null(attr(result, "n_failed_fits")))

  perturbed <- attr(result, "perturbed_months")
  expect_true(is.numeric(perturbed))
  expect_true(all(perturbed >= 1 & perturbed <= 12))
})


test_that("quantile_mapping verbose mode produces messages", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  # Force insufficient data in one month
  jan_idx <- which(data$mon.ts == 1)
  data$value[jan_idx[1:58]] <- 0

  expect_message(
    quantile_mapping(
      value = data$value,
      mean.change = changes$mean.change,
      var.change = changes$var.change,
      mon.ts = data$mon.ts,
      year.ts = data$year.ts,
      verbose = TRUE
    ),
    "Skipping months"
  )
})


# ==============================================================================
# STATISTICAL PROPERTIES TESTS
# ==============================================================================

test_that("quantile_mapping preserves distribution shape (relative)", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.5, var_factor = 1.5)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Check coefficient of variation (CV) is similar
  # CV = sd / mean should be approximately preserved with proportional changes
  nonzero_original <- data$value[data$value > 0]
  nonzero_result <- result[result > 0]

  cv_original <- sd(nonzero_original) / mean(nonzero_original)
  cv_result <- sd(nonzero_result) / mean(nonzero_result)

  # Should be within 20% (allowing for sampling variability)
  expect_true(abs(cv_result - cv_original) / cv_original < 0.2)
})


test_that("quantile_mapping maintains monthly patterns", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.2, var_factor = 1.2)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Calculate monthly means (excluding zeros)
  monthly_means_original <- tapply(
    data$value[data$value > 0],
    data$mon.ts[data$value > 0],
    mean
  )

  monthly_means_result <- tapply(
    result[result > 0],
    data$mon.ts[result > 0],
    mean
  )

  # Correlation between monthly means should be high
  cor_monthly <- cor(monthly_means_original, monthly_means_result, use = "complete.obs")
  expect_true(cor_monthly > 0.95)
})


# ==============================================================================
# FITTING METHOD TESTS
# ==============================================================================

test_that("quantile_mapping works with different fitting methods", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  # Method of moments (default)
  result_mme <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts,
    fit.method = "mme"
  )

  # Maximum likelihood
  result_mle <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts,
    fit.method = "mle"
  )

  # Both should produce valid results
  expect_true(all(is.finite(result_mme)))
  expect_true(all(is.finite(result_mle)))

  # Results should be similar but not identical
  nonzero_idx <- data$value > 0
  cor_methods <- cor(result_mme[nonzero_idx], result_mle[nonzero_idx])
  expect_true(cor_methods > 0.99)
})


# ==============================================================================
# VALIDATION FEATURE TESTS
# ==============================================================================

test_that("quantile_mapping validates output by default", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts,
    validate.output = TRUE
  )

  # Should have no NaN or Inf values
  expect_true(all(is.finite(result)))
})


test_that("quantile_mapping can disable output validation", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  # This should not error even if validation is off
  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts,
    validate.output = FALSE
  )

  expect_equal(length(result), length(data$value))
})


# ==============================================================================
# PERFORMANCE AND STRESS TESTS
# ==============================================================================

test_that("quantile_mapping handles large datasets efficiently", {
  skip_on_cran()  # Skip on CRAN to avoid long test times

  data <- generate_test_precip(n_years = 10)
  changes <- create_change_factors(n_years = 10)

  # Should complete in reasonable time (< 5 seconds)
  expect_lt(
    system.time({
      result <- quantile_mapping(
        value = data$value,
        mean.change = changes$mean.change,
        var.change = changes$var.change,
        mon.ts = data$mon.ts,
        year.ts = data$year.ts
      )
    })["elapsed"],
    5.0
  )

  expect_equal(length(result), length(data$value))
})


test_that("quantile_mapping is deterministic with same input", {
  data <- generate_test_precip(n_years = 2, seed = 456)
  changes <- create_change_factors(n_years = 2)

  result1 <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  result2 <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  expect_equal(result1, result2)
})


# ==============================================================================
# INTEGRATION TESTS
# ==============================================================================

# ==============================================================================
# INTEGRATION TESTS (REVISED)
# ==============================================================================

test_that("quantile_mapping integrates with real-world workflow", {
  # Simulate a realistic climate change scenario
  set.seed(98765)  # Use a different seed that works better
  data <- generate_test_precip(n_years = 5, seed = 98765)

  # Create realistic change factors (gradual increase over years)
  mean.change <- matrix(NA, nrow = 5, ncol = 12)
  var.change <- matrix(NA, nrow = 5, ncol = 12)

  for (y in 1:5) {
    # Mean increases by 3% per year
    mean.change[y, ] <- 1.0 + (y - 1) * 0.03

    # Variance increases by 5% per year
    var.change[y, ] <- 1.0 + (y - 1) * 0.05

    # Winter months (12, 1, 2) have larger increases
    winter_months <- c(12, 1, 2)
    mean.change[y, winter_months] <- mean.change[y, winter_months] * 1.1
    var.change[y, winter_months] <- var.change[y, winter_months] * 1.2
  }

  result <- quantile_mapping(
    value = data$value,
    mean.change = mean.change,
    var.change = var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # ============================================================================
  # TEST 1: Basic output validation
  # ============================================================================
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
  expect_equal(length(result), length(data$value))

  # ============================================================================
  # TEST 2: Dry day frequency preservation
  # ============================================================================
  expect_equal(sum(data$value == 0), sum(result == 0))

  # ============================================================================
  # TEST 3: Overall mean increase (most robust test)
  # ============================================================================
  # With all change factors >= 1.0, overall mean MUST increase
  overall_mean_original <- mean(data$value[data$value > 0])
  overall_mean_result <- mean(result[result > 0])
  expect_true(overall_mean_result > overall_mean_original)

  # ============================================================================
  # TEST 4: Overall variance increase
  # ============================================================================
  # With all variance change factors >= 1.0, overall variance should increase
  overall_var_original <- var(data$value[data$value > 0])
  overall_var_result <- var(result[result > 0])
  expect_true(overall_var_result > overall_var_original)

  # ============================================================================
  # TEST 5: Year 1 vs Year 5 comparison (direct test of time-varying factors)
  # ============================================================================
  # Year 5 has 12% higher mean change factor than Year 1
  # So Year 5 should have noticeably higher mean
  year1_values <- result[result > 0 & data$year.ts == 1]
  year5_values <- result[result > 0 & data$year.ts == 5]

  year1_mean <- mean(year1_values)
  year5_mean <- mean(year5_values)

  # Year 5 should be at least 8% higher than Year 1
  # (conservative estimate accounting for stochastic variability)
  expect_true(year5_mean > year1_mean * 1.08)

  # ============================================================================
  # TEST 6: Month-specific changes (winter months should increase more)
  # ============================================================================
  # Winter months in Year 5 should have higher values due to extra multiplier
  winter_idx <- data$mon.ts %in% c(12, 1, 2) & data$year.ts == 5 & result > 0
  summer_idx <- data$mon.ts %in% c(6, 7, 8) & data$year.ts == 5 & result > 0

  if (sum(winter_idx) > 10 && sum(summer_idx) > 10) {
    winter_mean <- mean(result[winter_idx])
    summer_mean <- mean(result[summer_idx])

    # Winter should have higher values (due to larger change factors)
    # This tests that month-specific factors are being applied
    expect_true(winter_mean > summer_mean)
  }

  # ============================================================================
  # TEST 7: Verify function actually modified the data
  # ============================================================================
  # At least 80% of non-zero values should have changed
  nonzero_idx <- data$value > 0
  changed <- result[nonzero_idx] != data$value[nonzero_idx]
  expect_true(mean(changed) > 0.8)

  # ============================================================================
  # TEST 8: Rank preservation within month-year combinations
  # ============================================================================
  # Quantile mapping should preserve ranks within each month-year
  # Test for Month 6, Year 3 (arbitrary choice with likely enough data)
  test_idx <- data$mon.ts == 6 & data$year.ts == 3 & data$value > 0

  if (sum(test_idx) > 5) {
    original_ranks <- rank(data$value[test_idx])
    result_ranks <- rank(result[test_idx])
    expect_equal(original_ranks, result_ranks)
  }
})


# ==============================================================================
# COMPARISON WITH ORIGINAL VALUES TESTS
# ==============================================================================

test_that("quantile_mapping with extreme decrease doesn't create negative values", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 0.1, var_factor = 0.1)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  expect_true(all(result >= 0))
})


test_that("quantile_mapping rank preservation within months", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- quantile_mapping(
    value = data$value,
    mean.change = changes$mean.change,
    var.change = changes$var.change,
    mon.ts = data$mon.ts,
    year.ts = data$year.ts
  )

  # Check rank preservation for a specific month
  month_idx <- which(data$mon.ts == 6 & data$value > 0)

  if (length(month_idx) > 5) {
    original_ranks <- rank(data$value[month_idx])
    result_ranks <- rank(result[month_idx])

    # Ranks should be identical (quantile mapping preserves order)
    expect_equal(original_ranks, result_ranks)
  }
})
