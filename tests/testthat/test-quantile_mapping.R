# ==============================================================================
# UNIT TESTS FOR perturb_prcp_qm()
# ==============================================================================

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
  prcp <- rgamma(n_days, shape = 2, scale = seasonal_mean / 2)

  # Add some dry days (~30%)
  dry_days <- sample(n_days, size = floor(n_days * 0.3))
  prcp[dry_days] <- 0

  list(
    prcp = prcp,
    month = month_idx,
    year = year_idx,
    n_days = n_days,
    n_years = n_years
  )
}


#' Create change factor matrices
#' @keywords internal
create_change_factors <- function(n_years, mean_factor = 1.1, var_factor = 1.2) {
  list(
    mean_factor = matrix(mean_factor, nrow = n_years, ncol = 12),
    var_factor  = matrix(var_factor,  nrow = n_years, ncol = 12)
  )
}


# ==============================================================================
# BASIC FUNCTIONALITY TESTS
# ==============================================================================

test_that("perturb_prcp_qm returns correct length output", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  expect_equal(length(result), length(data$prcp))
  expect_type(result, "double")
})


test_that("perturb_prcp_qm preserves zeros (dry days)", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  expect_equal(sum(data$prcp == 0), sum(result == 0))
  expect_true(all(result[data$prcp == 0] == 0))
})


test_that("perturb_prcp_qm increases mean when mean_factor > 1", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.5, var_factor = 1.0)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  nonzero_original <- data$prcp[data$prcp > 0]
  nonzero_result <- result[result > 0]

  expect_true(mean(nonzero_result) > mean(nonzero_original))
})


test_that("perturb_prcp_qm increases variance when var_factor > 1", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.0, var_factor = 2.0)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  nonzero_original <- data$prcp[data$prcp > 0]
  nonzero_result <- result[result > 0]

  expect_true(var(nonzero_result) > var(nonzero_original))
})


test_that("perturb_prcp_qm with no change returns similar distribution", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.0, var_factor = 1.0)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  nonzero_idx <- data$prcp > 0
  relative_diff <- abs(result[nonzero_idx] - data$prcp[nonzero_idx]) / data$prcp[nonzero_idx]

  expect_true(median(relative_diff) < 0.05)
})


# ==============================================================================
# INPUT VALIDATION TESTS
# ==============================================================================

test_that("perturb_prcp_qm rejects NULL inputs", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  expect_error(
    perturb_prcp_qm(
      prcp = NULL,
      mean_factor = changes$mean_factor,
      var_factor = changes$var_factor,
      month = data$month,
      year = data$year
    ),
    "'prcp' must not be NULL"
  )

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = NULL,
      var_factor = changes$var_factor,
      month = data$month,
      year = data$year
    ),
    "'mean_factor' must not be NULL"
  )
})


test_that("perturb_prcp_qm rejects mismatched lengths", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = changes$mean_factor,
      var_factor = changes$var_factor,
      month = data$month[1:100],  # wrong length
      year = data$year
    ),
    "'month' must have same length as 'prcp'"
  )
})


test_that("perturb_prcp_qm rejects negative precipitation", {
  data <- generate_test_precip(n_years = 2)
  data$prcp[1] <- -1
  changes <- create_change_factors(n_years = 2)

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = changes$mean_factor,
      var_factor = changes$var_factor,
      month = data$month,
      year = data$year
    ),
    "'prcp' must be non-negative"
  )
})


test_that("perturb_prcp_qm rejects invalid month indices", {
  data <- generate_test_precip(n_years = 2)
  data$month[1] <- 13
  changes <- create_change_factors(n_years = 2)

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = changes$mean_factor,
      var_factor = changes$var_factor,
      month = data$month,
      year = data$year
    ),
    "'month' must contain integers between 1 and 12"
  )
})


test_that("perturb_prcp_qm rejects non-matrix change factors", {
  data <- generate_test_precip(n_years = 2)

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = rep(1.1, 24),  # vector, not matrix
      var_factor = matrix(1.2, nrow = 2, ncol = 12),
      month = data$month,
      year = data$year
    ),
    "'mean_factor' must be a matrix with nrow = n_years, ncol = 12"
  )
})


test_that("perturb_prcp_qm rejects wrong matrix dimensions", {
  data <- generate_test_precip(n_years = 2)

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = matrix(1.1, nrow = 2, ncol = 6),  # wrong columns
      var_factor = matrix(1.2, nrow = 2, ncol = 12),
      month = data$month,
      year = data$year
    ),
    "'mean_factor' must have 12 columns \\(one per month\\)"
  )

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = matrix(1.1, nrow = 3, ncol = 12), # wrong rows (n_years = 2)
      var_factor = matrix(1.2, nrow = 3, ncol = 12),
      month = data$month,
      year = data$year
    ),
    "'mean_factor' must have 2 rows \\(one per year\\)"
  )
})




test_that("perturb_prcp_qm rejects negative change factors", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)
  changes$mean_factor[1, 1] <- -0.5

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = changes$mean_factor,
      var_factor = changes$var_factor,
      month = data$month,
      year = data$year
    ),
    "'mean_factor' must contain positive values"
  )
})


test_that("perturb_prcp_qm rejects zero change factors", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)
  changes$var_factor[1, 1] <- 0

  expect_error(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = changes$mean_factor,
      var_factor = changes$var_factor,
      month = data$month,
      year = data$year
    ),
    "'var_factor' must contain positive values"
  )
})


# ==============================================================================
# EDGE CASES AND ROBUSTNESS TESTS
# ==============================================================================

test_that("perturb_prcp_qm handles all-zero input", {
  data <- generate_test_precip(n_years = 2)
  data$prcp[] <- 0
  changes <- create_change_factors(n_years = 2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year,
    verbose = FALSE
  )

  expect_equal(result, data$prcp)
  expect_true(all(result == 0))
})


test_that("perturb_prcp_qm handles insufficient data in some months", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  jan_idx <- which(data$month == 1)
  data$prcp[jan_idx[1:55]] <- 0

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year,
    min_events = 10,
    verbose = FALSE
  )

  expect_equal(length(result), length(data$prcp))

  skipped <- attr(result, "skipped_months")
  expect_true(1 %in% skipped)
})


test_that("perturb_prcp_qm handles very small precipitation values", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  small_idx <- sample(which(data$prcp > 0), 50)
  data$prcp[small_idx] <- runif(50, 0.001, 0.01)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  expect_equal(length(result), length(data$prcp))
  expect_true(all(is.finite(result)))
})


test_that("perturb_prcp_qm handles extreme change factors", {
  data <- generate_test_precip(n_years = 2)

  changes_large <- create_change_factors(n_years = 2, mean_factor = 5.0, var_factor = 10.0)

  result_large <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes_large$mean_factor,
    var_factor = changes_large$var_factor,
    month = data$month,
    year = data$year
  )

  expect_true(all(is.finite(result_large)))
  expect_true(mean(result_large) > mean(data$prcp))

  changes_small <- create_change_factors(n_years = 2, mean_factor = 0.2, var_factor = 0.5)

  result_small <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes_small$mean_factor,
    var_factor = changes_small$var_factor,
    month = data$month,
    year = data$year
  )

  expect_true(all(is.finite(result_small)))
  expect_true(mean(result_small[result_small > 0]) < mean(data$prcp[data$prcp > 0]))
})


test_that("perturb_prcp_qm handles single year data", {
  data <- generate_test_precip(n_years = 1)
  changes <- create_change_factors(n_years = 1)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  expect_equal(length(result), length(data$prcp))
  expect_true(all(is.finite(result)))
})


test_that("perturb_prcp_qm handles varying change factors by year", {
  data <- generate_test_precip(n_years = 3)

  mean_factor <- matrix(NA_real_, nrow = 3, ncol = 12)
  var_factor  <- matrix(NA_real_, nrow = 3, ncol = 12)

  for (y in 1:3) {
    mean_factor[y, ] <- 1.0 + (y - 1) * 0.05
    var_factor[y, ]  <- 1.0 + (y - 1) * 0.10
  }

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = data$month,
    year = data$year
  )

  year_means <- tapply(result[result > 0], data$year[result > 0], mean)

  expect_true(year_means[2] > year_means[1])
  expect_true(year_means[3] > year_means[2])
})


# ==============================================================================
# ATTRIBUTE AND DIAGNOSTIC TESTS
# ==============================================================================

test_that("perturb_prcp_qm returns diagnostic attributes", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  expect_true(!is.null(attr(result, "perturbed_months")))
  expect_true(!is.null(attr(result, "skipped_months")))
  expect_true(!is.null(attr(result, "n_failed_fits")))

  perturbed <- attr(result, "perturbed_months")
  expect_true(is.numeric(perturbed))
  expect_true(all(perturbed >= 1 & perturbed <= 12))
})

test_that("perturb_prcp_qm diagnostics returns diagnostics object", {
  data <- generate_test_precip(n_years = 1)
  changes <- create_change_factors(n_years = 1)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year,
    diagnostics = TRUE
  )

  expect_type(result, "list")
  expect_true(all(c("adjusted", "diagnostics") %in% names(result)))
  expect_true(inherits(result$diagnostics, "prcp_qm_diagnostics"))
})


test_that("perturb_prcp_qm verbose mode produces messages", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  jan_idx <- which(data$month == 1)
  data$prcp[jan_idx[1:58]] <- 0

  expect_message(
    perturb_prcp_qm(
      prcp = data$prcp,
      mean_factor = changes$mean_factor,
      var_factor = changes$var_factor,
      month = data$month,
      year = data$year,
      verbose = TRUE
    ),
    "Skipping months"
  )
})


# ==============================================================================
# STATISTICAL PROPERTIES TESTS
# ==============================================================================

test_that("perturb_prcp_qm preserves distribution shape (relative)", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.5, var_factor = 1.5)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  nonzero_original <- data$prcp[data$prcp > 0]
  nonzero_result <- result[result > 0]

  cv_original <- sd(nonzero_original) / mean(nonzero_original)
  cv_result <- sd(nonzero_result) / mean(nonzero_result)

  expect_true(abs(cv_result - cv_original) / cv_original < 0.2)
})


test_that("perturb_prcp_qm maintains monthly patterns", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 1.2, var_factor = 1.2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  monthly_means_original <- tapply(
    data$prcp[data$prcp > 0],
    data$month[data$prcp > 0],
    mean
  )

  monthly_means_result <- tapply(
    result[result > 0],
    data$month[result > 0],
    mean
  )

  cor_monthly <- cor(monthly_means_original, monthly_means_result, use = "complete.obs")
  expect_true(cor_monthly > 0.95)
})


# ==============================================================================
# FITTING METHOD TESTS
# ==============================================================================

test_that("perturb_prcp_qm works with different fitting methods", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result_mme <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year,
    fit_method = "mme"
  )

  result_mle <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year,
    fit_method = "mle"
  )

  expect_true(all(is.finite(result_mme)))
  expect_true(all(is.finite(result_mle)))

  nonzero_idx <- data$prcp > 0
  cor_methods <- cor(result_mme[nonzero_idx], result_mle[nonzero_idx])
  expect_true(cor_methods > 0.99)
})


# ==============================================================================
# VALIDATION FEATURE TESTS
# ==============================================================================

test_that("perturb_prcp_qm validates output by default", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year,
    validate_output = TRUE
  )

  expect_true(all(is.finite(result)))
})


test_that("perturb_prcp_qm can disable output validation", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year,
    validate_output = FALSE
  )

  expect_equal(length(result), length(data$prcp))
})


# ==============================================================================
# PERFORMANCE AND STRESS TESTS
# ==============================================================================

test_that("perturb_prcp_qm handles large datasets efficiently", {
  skip_on_cran()

  data <- generate_test_precip(n_years = 10)
  changes <- create_change_factors(n_years = 10)

  expect_lt(
    system.time({
      result <- perturb_prcp_qm(
        prcp = data$prcp,
        mean_factor = changes$mean_factor,
        var_factor = changes$var_factor,
        month = data$month,
        year = data$year
      )
    })["elapsed"],
    5.0
  )

  expect_equal(length(result), length(data$prcp))
})


test_that("perturb_prcp_qm is deterministic with same input", {
  data <- generate_test_precip(n_years = 2, seed = 456)
  changes <- create_change_factors(n_years = 2)

  result1 <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  result2 <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  expect_equal(result1, result2)
})


# ==============================================================================
# INTEGRATION TESTS (REVISED)
# ==============================================================================

test_that("perturb_prcp_qm integrates with real-world workflow", {
  set.seed(98765)
  data <- generate_test_precip(n_years = 5, seed = 98765)

  mean_factor <- matrix(NA_real_, nrow = 5, ncol = 12)
  var_factor  <- matrix(NA_real_, nrow = 5, ncol = 12)

  for (y in 1:5) {
    mean_factor[y, ] <- 1.0 + (y - 1) * 0.03
    var_factor[y, ]  <- 1.0 + (y - 1) * 0.05

    winter_months <- c(12, 1, 2)
    mean_factor[y, winter_months] <- mean_factor[y, winter_months] * 1.1
    var_factor[y, winter_months]  <- var_factor[y, winter_months] * 1.2
  }

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = mean_factor,
    var_factor = var_factor,
    month = data$month,
    year = data$year
  )

  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
  expect_equal(length(result), length(data$prcp))

  expect_equal(sum(data$prcp == 0), sum(result == 0))

  overall_mean_original <- mean(data$prcp[data$prcp > 0])
  overall_mean_result <- mean(result[result > 0])
  expect_true(overall_mean_result > overall_mean_original)

  overall_var_original <- var(data$prcp[data$prcp > 0])
  overall_var_result <- var(result[result > 0])
  expect_true(overall_var_result > overall_var_original)

  year1_values <- result[result > 0 & data$year == 1]
  year5_values <- result[result > 0 & data$year == 5]

  year1_mean <- mean(year1_values)
  year5_mean <- mean(year5_values)

  expect_true(year5_mean > year1_mean * 1.08)

  winter_idx <- data$month %in% c(12, 1, 2) & data$year == 5 & result > 0
  summer_idx <- data$month %in% c(6, 7, 8) & data$year == 5 & result > 0

  if (sum(winter_idx) > 10 && sum(summer_idx) > 10) {
    winter_mean <- mean(result[winter_idx])
    summer_mean <- mean(result[summer_idx])
    expect_true(winter_mean > summer_mean)
  }

  nonzero_idx <- data$prcp > 0
  changed <- result[nonzero_idx] != data$prcp[nonzero_idx]
  expect_true(mean(changed) > 0.8)

  test_idx <- data$month == 6 & data$year == 3 & data$prcp > 0
  if (sum(test_idx) > 5) {
    original_ranks <- rank(data$prcp[test_idx])
    result_ranks <- rank(result[test_idx])
    expect_equal(original_ranks, result_ranks)
  }
})


# ==============================================================================
# COMPARISON WITH ORIGINAL VALUES TESTS
# ==============================================================================

test_that("perturb_prcp_qm with extreme decrease doesn't create negative values", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2, mean_factor = 0.1, var_factor = 0.1)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  expect_true(all(result >= 0))
})


test_that("perturb_prcp_qm rank preservation within months", {
  data <- generate_test_precip(n_years = 2)
  changes <- create_change_factors(n_years = 2)

  result <- perturb_prcp_qm(
    prcp = data$prcp,
    mean_factor = changes$mean_factor,
    var_factor = changes$var_factor,
    month = data$month,
    year = data$year
  )

  month_idx <- which(data$month == 6 & data$prcp > 0)

  if (length(month_idx) > 5) {
    original_ranks <- rank(data$prcp[month_idx])
    result_ranks <- rank(result[month_idx])
    expect_equal(original_ranks, result_ranks)
  }
})
