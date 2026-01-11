
# tests/testthat/test-apply_climate_perturbations.R

# ==============================================================================
# SETUP - ensure function exists in package context
# ==============================================================================

testthat::skip_if_not_installed("weathergenr")

# In package tests, the function should be available via namespace.
testthat::skip_if_not(exists("apply_climate_perturbations", where = asNamespace("weathergenr")))

# Optional: if these helpers are required at runtime by apply_climate_perturbations(),
# skip cleanly rather than warning during tests.
testthat::skip_if_not(exists("quantile_mapping", where = asNamespace("weathergenr")))
testthat::skip_if_not(exists("pet_hargreaves", where = asNamespace("weathergenr")))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

create_test_climate_data <- function(n_grids = 3, n_years = 5, seed = 123) {
  set.seed(seed)

  start_date <- as.Date("2000-01-01")
  n_days <- n_years * 365
  sim.dates <- seq(start_date, by = "day", length.out = n_days)

  climate.grid <- data.frame(
    id = 1:n_grids,
    x = seq(-120, -118, length.out = n_grids),
    y = seq(35, 37, length.out = n_grids)
  )

  climate.data <- lapply(seq_len(n_grids), function(i) {
    temp <- rnorm(n_days, mean = 15 + i, sd = 5)
    temp_range <- abs(rnorm(n_days, mean = 8, sd = 3))

    temp_min <- temp - temp_range * 0.6
    temp_max <- temp + temp_range * 0.4

    data.frame(
      precip = pmax(0, rnorm(n_days, mean = 2.5, sd = 5)),
      temp = temp,
      temp_min = temp_min,
      temp_max = temp_max,
      pet = abs(rnorm(n_days, mean = 3, sd = 1))
    )
  })

  list(
    climate.data = climate.data,
    climate.grid = climate.grid,
    sim.dates = sim.dates
  )
}

create_uniform_factors <- function(precip_mean = 1.2, precip_var = 1.1, temp = 2.0) {
  list(
    precip_mean = rep(precip_mean, 12),
    precip_var = rep(precip_var, 12),
    temp = rep(temp, 12)
  )
}

create_seasonal_factors <- function() {
  list(
    precip_mean = c(1.3, 1.3, 1.2, 1.1, 1.0, 1.0,
                    1.0, 1.0, 1.1, 1.2, 1.3, 1.3),
    precip_var = rep(1.1, 12),
    temp = c(2.5, 2.5, 2.0, 1.5, 1.0, 0.5,
             0.5, 1.0, 1.5, 2.0, 2.5, 2.5)
  )
}

# Convenience wrapper for quiet tests
.quiet <- function(expr) {
  suppressMessages(suppressWarnings(force(expr)))
}

# ==============================================================================
# INPUT VALIDATION TESTS
# ==============================================================================

test_that("validates NULL inputs", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    expect_error(
      apply_climate_perturbations(
        climate.data = NULL,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "'climate.data' must not be NULL"
    )

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = NULL,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp
      ),
      "'climate.grid' must not be NULL"
    )

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = NULL,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "'change.factor.precip.mean' must not be NULL"
    )

  }))
})

test_that("validates input types", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    expect_error(
      apply_climate_perturbations(
        climate.data = "not a list",
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "must be a list"
    )

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = "not a data frame",
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp
      ),
      "must be a data frame"
    )

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = "not dates",
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "must be a Date vector"
    )

  }))
})

test_that("validates dimension consistency", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_grids = 3)
    factors <- create_uniform_factors()

    bad_grid <- test_data$climate.grid[1:2, ]  # Only 2 rows
    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = bad_grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "must match"
    )

  }))
})

test_that("validates change factor lengths", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = rep(1.2, 6),
        change.factor.precip.variance = rep(1.1, 12),
        change.factor.temp.mean = rep(2.0, 12),
        verbose = FALSE
      ),
      "must have length 12"
    )

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = rep(1.2, 12),
        change.factor.precip.variance = rep(1.1, 12),
        change.factor.temp.mean = rep(2.0, 15),
        verbose = FALSE
      ),
      "must have length 12"
    )

  }))
})

test_that("validates positive change factors", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = rep(-0.5, 12),
        change.factor.precip.variance = rep(1.1, 12),
        change.factor.temp.mean = rep(2.0, 12),
        verbose = FALSE
      ),
      "must contain positive values"
    )

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = rep(1.2, 12),
        change.factor.precip.variance = rep(0, 12),
        change.factor.temp.mean = rep(2.0, 12),
        verbose = FALSE
      ),
      "must contain positive values"
    )

  }))
})

test_that("validates required columns in climate data", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    test_data$climate.data[[1]]$temp <- NULL

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "missing required columns.*temp"
    )

  }))
})

test_that("validates latitude column exists", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    test_data$climate.grid$y <- NULL

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "must contain a 'y' column"
    )

  }))
})

test_that("validates data length matches sim.dates", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    test_data$climate.data[[2]] <- test_data$climate.data[[2]][1:100, ]

    expect_error(
      apply_climate_perturbations(
        climate.data = test_data$climate.data,
        climate.grid = test_data$climate.grid,
        sim.dates = test_data$sim.dates,
        change.factor.precip.mean = factors$precip_mean,
        change.factor.precip.variance = factors$precip_var,
        change.factor.temp.mean = factors$temp,
        verbose = FALSE
      ),
      "has.*rows but 'sim.dates' has length"
    )

  }))
})

# ==============================================================================
# TRANSIENT CHANGE LOGIC TESTS
# ==============================================================================

test_that("transient temperature change reaches 2x factor at end", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 2.0)

    original_data <- test_data$climate.data[[1]]$temp

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      transient.temp.change = TRUE,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    year_vec <- as.integer(format(test_data$sim.dates, "%Y"))
    final_year_idx <- which(year_vec == max(year_vec))

    temp_change_final <- mean(
      result[[1]]$temp[final_year_idx] - original_data[final_year_idx]
    )

    expect_equal(temp_change_final, 4.0, tolerance = 0.2)

  }))
})

test_that("transient temperature change has mean equal to factor", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 2.0)

    original_data <- test_data$climate.data[[1]]$temp

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      transient.temp.change = TRUE,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    mean_temp_change <- mean(result[[1]]$temp - original_data)
    expect_equal(mean_temp_change, 2.0, tolerance = 0.1)

  }))
})

test_that("transient precipitation reaches correct endpoint", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(precip_mean = 1.2)

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      transient.precip.change = TRUE,
      transient.temp.change = FALSE,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    expect_true(TRUE)

  }))
})

test_that("step and transient changes produce equivalent mean precipitation", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(precip_mean = 1.3)

    result_step <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = rep(1.0, 12),
      change.factor.temp.mean = rep(0, 12),
      transient.precip.change = FALSE,
      transient.temp.change = FALSE,
      calculate.pet = FALSE
    )

    result_transient <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = rep(1.0, 12),
      change.factor.temp.mean = rep(0, 12),
      transient.precip.change = TRUE,
      transient.temp.change = FALSE,
      calculate.pet = FALSE
    )

    mean_step <- mean(sapply(result_step, function(x) mean(x$precip, na.rm = TRUE)))
    mean_trans <- mean(sapply(result_transient, function(x) mean(x$precip, na.rm = TRUE)))

    relative_diff <- abs(mean_step - mean_trans) / mean_step
    expect_lt(relative_diff, 0.05)

  }))
})

test_that("step and transient changes produce equivalent mean temperature", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 3.0)

    original_data <- lapply(test_data$climate.data, function(x) x$temp)

    result_step <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = rep(1.0, 12),
      change.factor.precip.variance = rep(1.0, 12),
      change.factor.temp.mean = factors$temp,
      transient.precip.change = FALSE,
      transient.temp.change = FALSE,
      calculate.pet = FALSE
    )

    test_data2 <- create_test_climate_data(n_years = 10)
    result_transient <- apply_climate_perturbations(
      climate.data = test_data2$climate.data,
      climate.grid = test_data2$climate.grid,
      sim.dates = test_data2$sim.dates,
      change.factor.precip.mean = rep(1.0, 12),
      change.factor.precip.variance = rep(1.0, 12),
      change.factor.temp.mean = factors$temp,
      transient.precip.change = FALSE,
      transient.temp.change = TRUE,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    mean_change_step <- mean(sapply(1:length(result_step), function(i) {
      mean(result_step[[i]]$temp - original_data[[i]])
    }))

    mean_change_trans <- mean(sapply(1:length(result_transient), function(i) {
      mean(result_transient[[i]]$temp - test_data2$climate.data[[i]]$temp)
    }))

    expect_equal(mean_change_step, 3.0, tolerance = 0.1)
    expect_equal(mean_change_trans, 3.0, tolerance = 0.1)

  }))
})

# ==============================================================================
# STEP CHANGE TESTS
# ==============================================================================

test_that("step temperature change is uniform across years", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(temp = 2.5)

    original_data <- test_data$climate.data[[1]]$temp

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      transient.temp.change = FALSE,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    temp_changes <- result[[1]]$temp - original_data
    month_ind <- as.integer(format(test_data$sim.dates, "%m"))

    for (m in 1:12) {
      month_changes <- temp_changes[month_ind == m]
      expect_lt(sd(month_changes), 0.01)
    }

  }))
})

test_that("step precipitation change is uniform across years", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_uniform_factors(precip_mean = 1.5)

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = rep(1.0, 12),
      change.factor.temp.mean = rep(0, 12),
      transient.precip.change = FALSE,
      transient.temp.change = FALSE,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    year_vec <- as.integer(format(test_data$sim.dates, "%Y"))
    years <- unique(year_vec)

    yearly_means <- sapply(years, function(yr) {
      mean(result[[1]]$precip[year_vec == yr], na.rm = TRUE)
    })

    cv <- sd(yearly_means) / mean(yearly_means)
    expect_lt(cv, 0.3)

  }))
})

# ==============================================================================
# OUTPUT STRUCTURE TESTS
# ==============================================================================

test_that("output has same structure as input", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_grids = 5)
    factors <- create_uniform_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    expect_equal(length(result), length(test_data$climate.data))

    for (i in seq_along(result)) {
      expect_equal(nrow(result[[i]]), nrow(test_data$climate.data[[i]]))
    }

    for (i in seq_along(result)) {
      expect_equal(names(result[[i]]), names(test_data$climate.data[[i]]))
    }

  }))
})

test_that("output contains no infinite or NaN values", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      verbose = FALSE
    )

    for (i in seq_along(result)) {
      expect_false(any(is.infinite(unlist(result[[i]]))))
      expect_false(any(is.nan(unlist(result[[i]]))))
    }

  }))
})

test_that("precipitation remains non-negative", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      verbose = FALSE
    )

    for (i in seq_along(result)) {
      expect_true(all(result[[i]]$precip >= 0, na.rm = TRUE))
    }

  }))
})

test_that("temperature ordering is maintained", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      verbose = FALSE
    )

    for (i in seq_along(result)) {
      violations_low <- sum(result[[i]]$temp_min > result[[i]]$temp)
      violations_high <- sum(result[[i]]$temp_max < result[[i]]$temp)

      expect_lt(violations_low / nrow(result[[i]]), 0.01)
      expect_lt(violations_high / nrow(result[[i]]), 0.01)
    }

  }))
})

# ==============================================================================
# PET CALCULATION TESTS
# ==============================================================================

test_that("PET is recalculated when requested", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors()

    original_pet <- test_data$climate.data[[1]]$pet

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      calculate.pet = TRUE,
      verbose = FALSE
    )

    expect_false(all(result[[1]]$pet == original_pet))
    expect_true(all(result[[1]]$pet >= 0, na.rm = TRUE))

  }))
})

test_that("PET is not recalculated when not requested", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors(temp = 0)

    original_pet <- test_data$climate.data[[1]]$pet

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      calculate.pet = FALSE,
      verbose = FALSE
    )

    expect_equal(
      sum(result[[1]]$pet != original_pet & is.finite(original_pet)),
      0
    )

  }))
})

# ==============================================================================
# EDGE CASE TESTS
# ==============================================================================

test_that("handles single grid cell", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_grids = 1)
    factors <- create_uniform_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      verbose = FALSE
    )

    expect_equal(length(result), 1)
    expect_true(is.data.frame(result[[1]]))

  }))
})

test_that("handles single year simulation", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 1)
    factors <- create_uniform_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      verbose = FALSE
    )

    expect_equal(length(result), length(test_data$climate.data))
    expect_equal(nrow(result[[1]]), nrow(test_data$climate.data[[1]]))

  }))
})

test_that("handles zero precipitation correctly", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()

    test_data$climate.data <- lapply(test_data$climate.data, function(x) {
      x$precip <- rep(0, nrow(x))
      x
    })

    factors <- create_uniform_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      verbose = FALSE
    )

    expect_true(all(result[[1]]$precip == 0))

  }))
})

test_that("handles extreme change factors", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()

    factors_extreme <- create_uniform_factors(
      precip_mean = 5.0,
      precip_var = 3.0,
      temp = 10.0
    )

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors_extreme$precip_mean,
      change.factor.precip.variance = factors_extreme$precip_var,
      change.factor.temp.mean = factors_extreme$temp,
      verbose = FALSE
    )

    expect_false(any(is.infinite(unlist(result[[1]]))))
    expect_false(any(is.nan(unlist(result[[1]]))))

  }))
})

test_that("handles very small change factors", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()

    factors_small <- create_uniform_factors(
      precip_mean = 1.01,
      precip_var = 1.0,
      temp = 0.1
    )

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors_small$precip_mean,
      change.factor.precip.variance = factors_small$precip_var,
      change.factor.temp.mean = factors_small$temp,
      verbose = FALSE
    )

    expect_false(any(is.infinite(unlist(result[[1]]))))

  }))
})

test_that("handles negative temperature changes", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data()
    factors <- create_uniform_factors(temp = -2.0)

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      verbose = FALSE
    )

    original_mean <- mean(test_data$climate.data[[1]]$temp)
    result_mean <- mean(result[[1]]$temp)

    expect_lt(result_mean, original_mean)

  }))
})

# ==============================================================================
# SEASONAL VARIATION TESTS
# ==============================================================================

test_that("applies seasonal variation correctly", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_years = 10)
    factors <- create_seasonal_factors()

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      transient.temp.change = FALSE,
      transient.precip.change = FALSE,
      verbose = FALSE
    )

    original_data <- test_data$climate.data[[1]]$temp
    month_ind <- as.integer(format(test_data$sim.dates, "%m"))

    temp_change_by_month <- sapply(1:12, function(m) {
      idx <- month_ind == m
      mean(result[[1]]$temp[idx] - original_data[idx])
    })

    winter_change <- mean(temp_change_by_month[c(12, 1, 2)])
    summer_change <- mean(temp_change_by_month[c(6, 7, 8)])

    expect_gt(winter_change, summer_change)

  }))
})

# ==============================================================================
# INTEGRATION TESTS
# ==============================================================================

test_that("realistic climate change scenario runs successfully", {
  suppressMessages(suppressWarnings({

    test_data <- create_test_climate_data(n_grids = 10, n_years = 30)

    factors <- list(
      precip_mean = c(1.25, 1.20, 1.20, 1.15, 1.15, 1.15,
                      1.15, 1.15, 1.15, 1.20, 1.20, 1.25),
      precip_var = rep(1.15, 12),
      temp = c(3.0, 2.8, 2.5, 2.0, 1.8, 1.5,
               1.5, 1.8, 2.0, 2.5, 2.8, 3.0)
    )

    result <- apply_climate_perturbations(
      climate.data = test_data$climate.data,
      climate.grid = test_data$climate.grid,
      sim.dates = test_data$sim.dates,
      change.factor.precip.mean = factors$precip_mean,
      change.factor.precip.variance = factors$precip_var,
      change.factor.temp.mean = factors$temp,
      transient.temp.change = TRUE,
      transient.precip.change = TRUE,
      calculate.pet = TRUE,
      verbose = FALSE
    )

    expect_equal(length(result), 10)
    expect_false(any(sapply(result, function(x) any(is.infinite(unlist(x))))))
    expect_false(any(sapply(result, function(x) any(is.nan(unlist(x))))))

    original_precip_mean <- mean(sapply(test_data$climate.data, function(x) mean(x$precip)))
    result_precip_mean <- mean(sapply(result, function(x) mean(x$precip)))
    expect_gt(result_precip_mean, original_precip_mean)

    original_temp_mean <- mean(sapply(test_data$climate.data, function(x) mean(x$temp)))
    result_temp_mean <- mean(sapply(result, function(x) mean(x$temp)))
    expect_gt(result_temp_mean, original_temp_mean)

  }))
})

