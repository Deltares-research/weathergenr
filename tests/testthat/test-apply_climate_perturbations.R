# tests/testthat/test-apply_climate_perturbations.R

# ==============================================================================
# SETUP - ensure function exists in package context
# ==============================================================================

testthat::skip_if_not_installed("weathergenr")

.ns <- asNamespace("weathergenr")
testthat::skip_if_not(exists("apply_climate_perturbations", where = .ns))
testthat::skip_if_not(exists("perturb_prcp_qm", where = .ns))
testthat::skip_if_not(exists("calculate_monthly_pet", where = .ns))

# Convenience: namespace-qualified call target
apply_cp <- weathergenr::apply_climate_perturbations

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

create_test_climate_data <- function(n_grids = 3, n_years = 5, seed = 123) {

  set.seed(seed)

  start_date <- as.Date("2000-01-01")
  n_days <- n_years * 365
  date <- seq(start_date, by = "day", length.out = n_days)

  grid <- data.frame(
    id  = 1:n_grids,
    x   = seq(-120, -118, length.out = n_grids),
    y   = seq(35, 37, length.out = n_grids),
    lat = seq(35, 37, length.out = n_grids)
  )

  data <- lapply(seq_len(n_grids), function(i) {

    temp <- rnorm(n_days, mean = 15 + i, sd = 5)
    temp_range <- abs(rnorm(n_days, mean = 8, sd = 3))

    temp_min <- temp - temp_range * 0.6
    temp_max <- temp + temp_range * 0.4

    data.frame(
      prcp     = pmax(0, rnorm(n_days, mean = 2.5, sd = 5)),
      temp     = temp,
      temp_min = temp_min,
      temp_max = temp_max,
      pet      = abs(rnorm(n_days, mean = 3, sd = 1))
    )
  })

  list(data = data, grid = grid, date = date)
}

create_uniform_factors <- function(precip_mean = 1.2, precip_var = 1.1, temp = 2.0) {
  list(
    precip_mean = rep(precip_mean, 12),
    precip_var  = rep(precip_var, 12),
    temp        = rep(temp, 12)
  )
}

create_seasonal_factors <- function() {
  list(
    precip_mean = c(1.3, 1.3, 1.2, 1.1, 1.0, 1.0,
                    1.0, 1.0, 1.1, 1.2, 1.3, 1.3),
    precip_var  = rep(1.1, 12),
    temp        = c(2.5, 2.5, 2.0, 1.5, 1.0, 0.5,
                    0.5, 1.0, 1.5, 2.0, 2.5, 2.5)
  )
}

.quiet <- function(expr) suppressMessages(suppressWarnings(force(expr)))

# ==============================================================================
# INPUT VALIDATION TESTS
# ==============================================================================

testthat::test_that("validates NULL inputs", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    testthat::expect_error(
      apply_cp(
        data = NULL,
        grid = test_data$grid,
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        verbose = FALSE
      ),
      "climate\\.data.*must not be NULL"
    )

    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = NULL,
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        verbose = FALSE
      ),
      "'grid' must not be NULL"
    )

    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = test_data$grid,
        date = test_data$date,
        prcp_mean_factor = NULL,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        verbose = FALSE
      ),
      "change\\.factor\\.precip\\.mean.*must not be NULL"
    )

  })
})

testthat::test_that("validates input types", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    testthat::expect_error(
      apply_cp(
        data = "not a list",
        grid = test_data$grid,
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        verbose = FALSE
      ),
      "must be a list"
    )

    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = "not a data frame",
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        verbose = FALSE
      ),
      "must be a data frame"
    )

    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = test_data$grid,
        date = "not dates",
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        verbose = FALSE
      ),
      "must be a Date vector"
    )

  })
})

testthat::test_that("validates dimension consistency", {
  .quiet({

    test_data <- create_test_climate_data(n_grids = 3)
    factors <- create_uniform_factors()

    bad_grid <- test_data$grid[1:2, ]
    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = bad_grid,
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        verbose = FALSE
      ),
      "must match"
    )

  })
})

testthat::test_that("validates required columns in climate data", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    test_data$data[[1]]$temp <- NULL

    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = test_data$grid,
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        compute_pet      = FALSE,
        verbose = FALSE
      ),
      "missing columns:.*temp"
    )

  })
})

testthat::test_that("validates latitude column exists (lat)", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    test_data$grid$lat <- NULL

    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = test_data$grid,
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        compute_pet      = TRUE,
        pet_method       = "hargreaves",
        verbose = FALSE
      ),
      "lat|y"
    )

  })
})

testthat::test_that("validates data length matches date", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    test_data$data[[2]] <- test_data$data[[2]][1:100, ]

    testthat::expect_error(
      apply_cp(
        data = test_data$data,
        grid = test_data$grid,
        date = test_data$date,
        prcp_mean_factor = factors$precip_mean,
        prcp_var_factor  = factors$precip_var,
        temp_delta       = factors$temp,
        compute_pet      = FALSE,
        verbose = FALSE
      ),
      "row count does not match|sim\\.dates"
    )

  })
})

# ==============================================================================
# TRANSIENT CHANGE LOGIC TESTS
# ==============================================================================

testthat::test_that("transient temperature change reaches ~2x delta at end (mean equals delta)", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 2.0)

    original_data <- test_data$data[[1]]$temp

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      temp_transient   = TRUE,
      prcp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    year_vec <- as.integer(format(test_data$date, "%Y"))
    final_year_idx <- which(year_vec == max(year_vec))
    temp_change_final <- mean(result[[1]]$temp[final_year_idx] - original_data[final_year_idx])

    testthat::expect_equal(temp_change_final, 4.0, tolerance = 0.25)

  })
})

testthat::test_that("transient temperature change has mean approximately equal to delta", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 2.0)

    original_data <- test_data$data[[1]]$temp

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      temp_transient   = TRUE,
      prcp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    mean_temp_change <- mean(result[[1]]$temp - original_data)
    testthat::expect_equal(mean_temp_change, 2.0, tolerance = 0.15)

  })
})

testthat::test_that("step and transient changes produce similar mean precipitation", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(precip_mean = 1.3)

    result_step <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = rep(0, 12),
      prcp_transient   = FALSE,
      temp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    result_transient <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = rep(0, 12),
      prcp_transient   = TRUE,
      temp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    mean_step <- mean(vapply(result_step, function(x) mean(x$prcp, na.rm = TRUE), numeric(1)))
    mean_trans <- mean(vapply(result_transient, function(x) mean(x$prcp, na.rm = TRUE), numeric(1)))

    relative_diff <- abs(mean_step - mean_trans) / mean_step
    testthat::expect_lt(relative_diff, 0.05)

  })
})

testthat::test_that("step and transient temperature changes produce expected mean deltas", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 3.0)

    original_data <- lapply(test_data$data, function(x) x$temp)

    result_step <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      prcp_transient   = FALSE,
      temp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    test_data2 <- create_test_climate_data(n_years = 10)
    result_transient <- apply_cp(
      data = test_data2$data,
      grid = test_data2$grid,
      date = test_data2$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      prcp_transient   = FALSE,
      temp_transient   = TRUE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    mean_change_step <- mean(vapply(seq_along(result_step), function(i) {
      mean(result_step[[i]]$temp - original_data[[i]])
    }, numeric(1)))

    mean_change_trans <- mean(vapply(seq_along(result_transient), function(i) {
      mean(result_transient[[i]]$temp - test_data2$data[[i]]$temp)
    }, numeric(1)))

    testthat::expect_equal(mean_change_step, 3.0, tolerance = 0.15)
    testthat::expect_equal(mean_change_trans, 3.0, tolerance = 0.20)

  })
})

# ==============================================================================
# STEP CHANGE TESTS
# ==============================================================================

testthat::test_that("step temperature change is uniform within each month", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 2.5)

    original_data <- test_data$data[[1]]$temp

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      temp_transient   = FALSE,
      prcp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    temp_changes <- result[[1]]$temp - original_data
    month_ind <- as.integer(format(test_data$date, "%m"))

    for (m in 1:12) {
      month_changes <- temp_changes[month_ind == m]
      testthat::expect_lt(stats::sd(month_changes), 0.01)
    }

  })
})

testthat::test_that("step precipitation change is reasonably stable across years", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(precip_mean = 1.5)

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = rep(0, 12),
      prcp_transient   = FALSE,
      temp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    year_vec <- as.integer(format(test_data$date, "%Y"))
    years <- unique(year_vec)

    yearly_means <- vapply(years, function(yr) {
      mean(result[[1]]$prcp[year_vec == yr], na.rm = TRUE)
    }, numeric(1))

    cv <- stats::sd(yearly_means) / mean(yearly_means)
    testthat::expect_lt(cv, 0.5)

  })
})

# ==============================================================================
# OUTPUT STRUCTURE TESTS
# ==============================================================================

testthat::test_that("output has same structure as input", {
  .quiet({

    test_data <- create_test_climate_data(n_grids = 5)
    factors <- create_uniform_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    testthat::expect_equal(length(result), length(test_data$data))

    for (i in seq_along(result)) {
      testthat::expect_equal(nrow(result[[i]]), nrow(test_data$data[[i]]))
      testthat::expect_equal(names(result[[i]]), names(test_data$data[[i]]))
    }

  })
})

testthat::test_that("output contains no infinite or NaN values", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    for (i in seq_along(result)) {
      testthat::expect_false(any(is.infinite(unlist(result[[i]]))))
      testthat::expect_false(any(is.nan(unlist(result[[i]]))))
    }

  })
})

testthat::test_that("precipitation remains non-negative", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    for (i in seq_along(result)) {
      testthat::expect_true(all(result[[i]]$prcp >= 0, na.rm = TRUE))
    }

  })
})

testthat::test_that("temperature ordering is maintained", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    for (i in seq_along(result)) {
      violations_low  <- sum(result[[i]]$temp_min > result[[i]]$temp, na.rm = TRUE)
      violations_high <- sum(result[[i]]$temp_max < result[[i]]$temp, na.rm = TRUE)

      testthat::expect_lt(violations_low / nrow(result[[i]]), 0.01)
      testthat::expect_lt(violations_high / nrow(result[[i]]), 0.01)
    }

  })
})

# ==============================================================================
# PET CALCULATION TESTS
# ==============================================================================

testthat::test_that("PET is recalculated when requested", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors(temp = 1.0)

    original_pet <- test_data$data[[1]]$pet

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      compute_pet      = TRUE,
      pet_method       = "hargreaves",
      verbose = FALSE
    )

    testthat::expect_false(isTRUE(all.equal(result[[1]]$pet, original_pet)))
    testthat::expect_true(all(result[[1]]$pet >= 0, na.rm = TRUE))

  })
})

testthat::test_that("PET is not recalculated when not requested", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors(temp = 0)

    original_pet <- test_data$data[[1]]$pet

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    testthat::expect_equal(result[[1]]$pet, original_pet)

  })
})

# ==============================================================================
# EDGE CASE TESTS
# ==============================================================================

testthat::test_that("handles single grid cell", {
  .quiet({

    test_data <- create_test_climate_data(n_grids = 1)
    factors <- create_uniform_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    testthat::expect_equal(length(result), 1)
    testthat::expect_true(is.data.frame(result[[1]]))

  })
})

testthat::test_that("handles single year simulation", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 1)
    factors <- create_uniform_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    testthat::expect_equal(length(result), length(test_data$data))
    testthat::expect_equal(nrow(result[[1]]), nrow(test_data$data[[1]]))

  })
})

testthat::test_that("handles zero precipitation correctly", {
  .quiet({

    test_data <- create_test_climate_data()

    test_data$data <- lapply(test_data$data, function(x) {
      x$prcp <- rep(0, nrow(x))
      x
    })

    factors <- create_uniform_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    testthat::expect_true(all(result[[1]]$prcp == 0))

  })
})

testthat::test_that("handles negative temperature changes", {
  .quiet({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors(temp = -2.0)

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    original_mean <- mean(test_data$data[[1]]$temp)
    result_mean <- mean(result[[1]]$temp)
    testthat::expect_lt(result_mean, original_mean)

  })
})

# ==============================================================================
# SEASONAL VARIATION TESTS
# ==============================================================================

testthat::test_that("applies seasonal temperature variation correctly", {
  .quiet({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_seasonal_factors()

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = rep(1.0, 12),
      prcp_var_factor  = rep(1.0, 12),
      temp_delta       = factors$temp,
      temp_transient   = FALSE,
      prcp_transient   = FALSE,
      compute_pet      = FALSE,
      verbose = FALSE
    )

    original_data <- test_data$data[[1]]$temp
    month_ind <- as.integer(format(test_data$date, "%m"))

    temp_change_by_month <- vapply(1:12, function(m) {
      idx <- month_ind == m
      mean(result[[1]]$temp[idx] - original_data[idx])
    }, numeric(1))

    winter_change <- mean(temp_change_by_month[c(12, 1, 2)])
    summer_change <- mean(temp_change_by_month[c(6, 7, 8)])

    testthat::expect_gt(winter_change, summer_change)

  })
})

# ==============================================================================
# INTEGRATION TESTS
# ==============================================================================

testthat::test_that("realistic climate change scenario runs successfully", {
  .quiet({

    test_data <- create_test_climate_data(n_grids = 10, n_years = 30)

    factors <- list(
      precip_mean = c(1.25, 1.20, 1.20, 1.15, 1.15, 1.15,
                      1.15, 1.15, 1.15, 1.20, 1.20, 1.25),
      precip_var  = rep(1.15, 12),
      temp        = c(3.0, 2.8, 2.5, 2.0, 1.8, 1.5,
                      1.5, 1.8, 2.0, 2.5, 2.8, 3.0)
    )

    result <- apply_cp(
      data = test_data$data,
      grid = test_data$grid,
      date = test_data$date,
      prcp_mean_factor = factors$precip_mean,
      prcp_var_factor  = factors$precip_var,
      temp_delta       = factors$temp,
      temp_transient   = TRUE,
      prcp_transient   = TRUE,
      compute_pet      = TRUE,
      pet_method       = "hargreaves",
      verbose = FALSE
    )

    testthat::expect_equal(length(result), 10)
    testthat::expect_false(any(vapply(result, function(x) any(is.infinite(unlist(x))), logical(1))))
    testthat::expect_false(any(vapply(result, function(x) any(is.nan(unlist(x))),      logical(1))))

    original_prcp_mean <- mean(vapply(test_data$data, function(x) mean(x$prcp), numeric(1)))
    result_prcp_mean   <- mean(vapply(result,         function(x) mean(x$prcp), numeric(1)))
    testthat::expect_gt(result_prcp_mean, original_prcp_mean)

    original_temp_mean <- mean(vapply(test_data$data, function(x) mean(x$temp), numeric(1)))
    result_temp_mean   <- mean(vapply(result,         function(x) mean(x$temp), numeric(1)))
    testthat::expect_gt(result_temp_mean, original_temp_mean)

  })
})

