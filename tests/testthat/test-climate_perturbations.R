# ==============================================================================
# Unit Tests for climate_perturbations.R
# ==============================================================================
library(testthat)

# ------------------------------------------------------------------------------
# Helpers for deterministic synthetic inputs
# ------------------------------------------------------------------------------

.make_dates_noleap_years <- function(start = "2000-01-01", n_years = 2L) {
  as.Date(start) + 0:(365L * n_years - 1L)
}

# Build precip with:
# - Enough wet days per (year,month) to satisfy min_events=10
# - Non-degenerate wet-day distribution (Gamma-like) so fits don't collapse
.make_precip_with_min_wet <- function(date, min_wet = 15L, seed = 1L) {
  set.seed(seed)
  n <- length(date)
  month <- as.integer(format(date, "%m"))
  year  <- as.integer(format(date, "%Y"))
  year_idx <- year - min(year) + 1L

  precip <- numeric(n)

  keys <- unique(data.frame(y = year_idx, m = month))
  for (k in seq_len(nrow(keys))) {
    yk <- keys$y[k]
    mk <- keys$m[k]
    idx <- which(year_idx == yk & month == mk)
    n_idx <- length(idx)

    # Pick exactly min_wet wet days (or all days if month shorter in edge cases)
    n_wet <- min(min_wet, n_idx)
    wet_local <- sample(idx, n_wet, replace = FALSE)

    # Wet amounts with variance (avoid degenerate fit)
    precip[wet_local] <- rgamma(length(wet_local), shape = 2.0, scale = 3.0)

    # Remaining days stay dry (0)
  }

  precip
}

# ==============================================================================
# Core behavior tests
# ==============================================================================

test_that("apply_climate_perturbations: returns list-of-data.frames with same structure", {
  n_years <- 2L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1:2, lat = c(0, 10))

  mk_cell <- function(seed) {
    set.seed(seed)
    precip <- rgamma(n_days, shape = 1.2, scale = 5)
    precip[sample.int(n_days, 150)] <- 0
    temp <- 20 + sin(2 * pi * (1:n_days) / 365) * 3
    data.frame(
      precip = precip,
      temp = temp,
      temp_min = temp - 2,
      temp_max = temp + 2
    )
  }

  data <- list(mk_cell(1), mk_cell(2))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1, 12),
    precip_var_factor = rep(1, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    precip_transient = FALSE,
    precip_occurrence_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    diagnostic = FALSE,          # KEY: return list of data.frames
    seed = 1,
    verbose = FALSE
  )

  expect_true(is.list(out))
  expect_length(out, length(data))

  for (i in seq_along(out)) {
    expect_true(is.data.frame(out[[i]]))
    expect_equal(nrow(out[[i]]), n_days)
    expect_true(all(c("precip", "temp", "temp_min", "temp_max") %in% names(out[[i]])))
  }
})

test_that("apply_climate_perturbations: diagnostic=TRUE returns list with data and diagnostic", {
  n_years <- 1L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 0)

  data <- list(data.frame(
    precip = rgamma(n_days, shape = 1.5, scale = 4),
    temp = rep(10, n_days),
    temp_min = rep(9, n_days),
    temp_max = rep(11, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1, 12),
    precip_var_factor = rep(1, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    precip_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    diagnostic = TRUE,
    verbose = FALSE
  )

  expect_type(out, "list")
  expect_named(out, c("data", "diagnostic"))
  expect_true(is.list(out$data))
  expect_true(is.list(out$diagnostic))
  expect_true(is.data.frame(out$data[[1]]))
})

test_that("apply_climate_perturbations: temperature deltas apply (step change)", {
  n_years <- 2L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 0)

  precip <- rep(1, n_days)
  temp <- rep(10, n_days)

  data <- list(data.frame(
    precip = precip,
    temp = temp,
    temp_min = temp - 1,
    temp_max = temp + 1
  ))

  delta <- rep(2, 12)

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1, 12),
    precip_var_factor = rep(1, 12),
    temp_delta = delta,
    temp_transient = FALSE,
    precip_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    diagnostic = FALSE,
    verbose = FALSE
  )[[1]]

  expect_equal(out$temp, temp + 2)
  expect_equal(out$temp_min, (temp - 1) + 2)
  expect_equal(out$temp_max, (temp + 1) + 2)
})

test_that("apply_climate_perturbations: transient temperature ramp has 0 at first year and ~2*delta at last year", {
  n_years <- 3L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 0)

  temp <- rep(0, n_days)
  data <- list(data.frame(
    precip = rep(1, n_days),
    temp = temp,
    temp_min = temp,
    temp_max = temp
  ))

  delta <- rep(1, 12)

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1, 12),
    precip_var_factor = rep(1, 12),
    temp_delta = delta,
    temp_transient = TRUE,
    precip_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    diagnostic = FALSE,
    verbose = FALSE
  )[[1]]

  cal_year <- as.integer(format(date, "%Y"))
  year_idx <- cal_year - min(cal_year) + 1L

  i1 <- which(year_idx == 1 & format(date, "%m") == "01")[1]
  i3 <- which(year_idx == 3 & format(date, "%m") == "01")[1]

  expect_equal(out$temp[i1], 0, tolerance = 1e-12)
  expect_equal(out$temp[i3], 2, tolerance = 1e-12)
})

test_that("apply_climate_perturbations: precipitation unchanged under identity factors when enforce_target_mean=FALSE", {
  n_years <- 2L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 0)

  set.seed(3)
  precip <- rgamma(n_days, shape = 1.5, scale = 4)
  precip[sample.int(n_days, 120)] <- 0

  data <- list(data.frame(
    precip = precip,
    temp = rep(0, n_days),
    temp_min = rep(0, n_days),
    temp_max = rep(0, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1, 12),
    precip_var_factor = rep(1, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    precip_transient = FALSE,
    precip_occurrence_factor = NULL,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    diagnostic = FALSE,
    seed = 1,
    verbose = FALSE
  )[[1]]

  wet <- precip > 0
  expect_identical(out$precip[!wet], precip[!wet])
  expect_lt(mean(abs(out$precip[wet] - precip[wet])), 1e-6)
})

test_that("apply_climate_perturbations: occurrence factor increases wet-day count when Gamma fit is feasible", {
  n_years <- 2L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 0)

  precip <- .make_precip_with_min_wet(date, min_wet = 15L, seed = 7L)

  data <- list(data.frame(
    precip = precip,
    temp = rep(0, n_days),
    temp_min = rep(0, n_days),
    temp_max = rep(0, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1, 12),
    precip_var_factor = rep(1, 12),
    precip_occurrence_factor = rep(1.5, 12),
    precip_intensity_threshold = 0,
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    precip_transient = FALSE,
    precip_occurrence_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    diagnostic = TRUE,      # ensures diagnostics requested; occurrence also requests diagnostics internally
    seed = 42,
    verbose = FALSE
  )

  # unwrap
  out1 <- out$data[[1]]

  wet0 <- precip > 0
  wet1 <- out1$precip > 0

  expect_gte(sum(wet1), sum(wet0))

  added <- (!wet0) & wet1
  if (any(added)) {
    expect_true(all(out1$precip[added] > 0))
  }
})

test_that("apply_climate_perturbations: safety rails enforce cap and floor", {
  n_years <- 2L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 0)

  set.seed(9)
  precip <- rgamma(n_days, shape = 0.8, scale = 20)
  precip[sample.int(n_days, 100)] <- 0

  data <- list(data.frame(
    precip = precip,
    temp = rep(0, n_days),
    temp_min = rep(0, n_days),
    temp_max = rep(0, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1.2, 12),
    precip_var_factor = rep(1.0, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    precip_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = TRUE,
    extreme_k = 1.5,
    enforce_target_mean = TRUE,
    precip_floor_mm_day = 0.2,
    precip_cap_mm_day = 50,
    diagnostic = FALSE,
    seed = 1,
    verbose = FALSE
  )[[1]]

  wet <- out$precip > 0
  if (any(wet)) {
    expect_true(all(out$precip[wet] >= 0.2))
    expect_true(max(out$precip, na.rm = TRUE) <= 50)
  }
})

test_that("apply_climate_perturbations: PET is added when compute_pet=TRUE", {
  skip_if_not(exists("calculate_monthly_pet", where = asNamespace("weathergenr"), inherits = FALSE))

  n_years <- 2L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 10)

  data <- list(data.frame(
    precip = rep(1, n_days),
    temp = rep(20, n_days),
    temp_min = rep(18, n_days),
    temp_max = rep(22, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    precip_mean_factor = rep(1, 12),
    precip_var_factor = rep(1, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    precip_transient = FALSE,
    compute_pet = TRUE,
    pet_method = "hargreaves",
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    diagnostic = FALSE,
    verbose = FALSE
  )[[1]]

  expect_true("pet" %in% names(out))
  expect_length(out$pet, n_days)
  expect_true(all(is.finite(out$pet)))
  expect_true(all(out$pet >= 0))
})

# ==============================================================================
# Additional coverage (validation + scale_var_with_mean branch)
# ==============================================================================

test_that("apply_climate_perturbations: input validation errors are thrown", {
  date <- .make_dates_noleap_years("2000-01-01", 1L)
  grid <- data.frame(id = 1, lat = 0)
  data <- list(data.frame(
    precip = rep(1, length(date)),
    temp = rep(0, length(date)),
    temp_min = rep(0, length(date)),
    temp_max = rep(0, length(date))
  ))

  expect_error(apply_climate_perturbations(data = NULL, grid = grid, date = date,
                                           precip_mean_factor = rep(1, 12),
                                           precip_var_factor = rep(1, 12),
                                           temp_delta = rep(0, 12)),
               "'climate.data' must not be NULL")

  expect_error(apply_climate_perturbations(data = data, grid = NULL, date = date,
                                           precip_mean_factor = rep(1, 12),
                                           precip_var_factor = rep(1, 12),
                                           temp_delta = rep(0, 12)),
               "'grid' must not be NULL")

  expect_error(apply_climate_perturbations(data = data, grid = grid, date = NULL,
                                           precip_mean_factor = rep(1, 12),
                                           precip_var_factor = rep(1, 12),
                                           temp_delta = rep(0, 12)),
               "'sim.dates' must not be NULL")

  expect_error(apply_climate_perturbations(data = data, grid = grid, date = date,
                                           precip_mean_factor = rep(1, 11),
                                           precip_var_factor = rep(1, 12),
                                           temp_delta = rep(0, 12)),
               "precip_mean_factor")
})

test_that("apply_climate_perturbations: scale_var_with_mean ignores precip_var_factor", {
  n_years <- 1L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n_days <- length(date)

  grid <- data.frame(id = 1, lat = 0)

  data <- list(data.frame(
    precip = rgamma(n_days, shape = 1.5, scale = 3),
    temp = rep(0, n_days),
    temp_min = rep(0, n_days),
    temp_max = rep(0, n_days)
  ))

  expect_warning(
    apply_climate_perturbations(
      data = data,
      grid = grid,
      date = date,
      precip_mean_factor = rep(1.1, 12),
      precip_var_factor = rep(9, 12),     # should be ignored
      temp_delta = rep(0, 12),
      temp_transient = FALSE,
      precip_transient = FALSE,
      compute_pet = FALSE,
      scale_var_with_mean = TRUE,
      exaggerate_extremes = FALSE,
      enforce_target_mean = FALSE,
      diagnostic = FALSE,
      verbose = TRUE
    ),
    "Ignoring 'precip_var_factor'"
  )
})
