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
.handoff_cell <- function(date, seed = 5L) {
  n <- length(date)
  set.seed(seed)
  cell <- data.frame(
    precip   = .make_precip_with_min_wet(date, min_wet = 15L, seed = seed),
    temp     = stats::rnorm(n, 25, 3),
    temp_min = stats::rnorm(n, 20, 3)
  )
  cell$temp_max <- pmax(cell$temp_min + 5, stats::rnorm(n, 30, 3))
  cell
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

  # Strictly greater, not >=. The >= form passed while occurrence perturbation
  # was a no-op: apply_climate_perturbations() read base_gamma/target_gamma off
  # the quantile-mapping result, which never returned them, so the branch was
  # unreachable and wet1 always equalled wet0.
  expect_gt(sum(wet1), sum(wet0))

  # A 1.5 factor targets ~1.5x the wet days, subject to rounding per
  # (month, year) group and to running out of dry days to convert.
  expect_equal(sum(wet1) / sum(wet0), 1.5, tolerance = 0.1)

  added <- (!wet0) & wet1
  expect_true(any(added))
  expect_true(all(out1$precip[added] > 0))
})

test_that("occurrence factors below 1 remove wet days", {
  n_years <- 3L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  grid <- data.frame(id = 1L, lat = 0)
  dat <- list(.handoff_cell(date, seed = 9L))

  run <- function(f) {
    out <- apply_climate_perturbations(
      data = dat, grid = grid, date = date,
      precip_mean_factor = rep(1, 12), precip_var_factor = rep(1, 12),
      precip_occurrence_factor = rep(f, 12),
      temp_delta = rep(0, 12), temp_transient = FALSE,
      precip_transient = FALSE, precip_occurrence_transient = FALSE,
      scale_var_with_mean = FALSE, enforce_target_mean = FALSE,
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE, seed = 42L
    )
    mean(out[[1]]$precip > 0)
  }

  base_frac <- mean(dat[[1]]$precip > 0)

  expect_equal(run(1.0) / base_frac, 1.0, tolerance = 0.02)
  expect_equal(run(0.6) / base_frac, 0.6, tolerance = 0.05)
  expect_lt(run(0.6), run(1.0))
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
               "'data' must not be NULL")

  expect_error(apply_climate_perturbations(data = data, grid = NULL, date = date,
                                           precip_mean_factor = rep(1, 12),
                                           precip_var_factor = rep(1, 12),
                                           temp_delta = rep(0, 12)),
               "'grid' must not be NULL")

  expect_error(apply_climate_perturbations(data = data, grid = grid, date = NULL,
                                           precip_mean_factor = rep(1, 12),
                                           precip_var_factor = rep(1, 12),
                                           temp_delta = rep(0, 12)),
               "'date' must not be NULL")

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

# ==============================================================================
# apply_climate_perturbations(): remaining validation branches
#
# The existing validation test above covers the three NULL guards and the
# factor-length check. These cover the rest -- type, conformity, safety-rail
# and change-matrix checks -- all of which exit before any perturbation work,
# so the whole block costs milliseconds.
# ==============================================================================

.cp_valid_inputs <- function(n_years = 1L) {
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  n <- length(date)

  list(
    date = date,
    grid = data.frame(id = 1, lat = 0),
    data = list(data.frame(
      precip   = rep(1, n),
      temp     = rep(10, n),
      temp_min = rep(8, n),
      temp_max = rep(12, n)
    )),
    precip_mean_factor = rep(1, 12),
    precip_var_factor  = rep(1, 12),
    temp_delta         = rep(0, 12)
  )
}

.cp_call <- function(inputs, ...) {
  args <- list(
    data               = inputs$data,
    grid               = inputs$grid,
    date               = inputs$date,
    precip_mean_factor = inputs$precip_mean_factor,
    precip_var_factor  = inputs$precip_var_factor,
    temp_delta         = inputs$temp_delta
  )

  # Replace wholesale: modifyList() would merge the list-valued arguments.
  overrides <- list(...)
  args[names(overrides)] <- overrides

  do.call(apply_climate_perturbations, args)
}

test_that("apply_climate_perturbations: required change factors must be supplied", {
  inputs <- .cp_valid_inputs()

  expect_error(
    .cp_call(inputs, precip_mean_factor = NULL),
    "'precip_mean_factor' must not be NULL"
  )
  expect_error(
    .cp_call(inputs, temp_delta = NULL),
    "'temp_delta' must not be NULL"
  )

  # the variance factor is only required when it is not derived from the mean
  expect_error(
    .cp_call(inputs, precip_var_factor = NULL, scale_var_with_mean = FALSE),
    "'precip_var_factor' must not be NULL"
  )
})

test_that("apply_climate_perturbations: core argument types are checked", {
  inputs <- .cp_valid_inputs()

  expect_error(.cp_call(inputs, data = "not-a-list"), "'data' must be a list of data frames")
  expect_error(.cp_call(inputs, grid = list(lat = 0)), "'grid' must be a data frame")
  expect_error(
    .cp_call(inputs, date = as.character(inputs$date)),
    "'date' must be a Date vector"
  )
})

test_that("apply_climate_perturbations: data, grid and dates must conform", {
  inputs <- .cp_valid_inputs()

  # two grid rows but only one data cell
  expect_error(
    .cp_call(inputs, grid = data.frame(id = 1:2, lat = c(0, 10))),
    "must match.*number of rows in 'grid'"
  )

  # neither 'lat' nor 'y' present
  expect_error(
    .cp_call(inputs, grid = data.frame(id = 1, elevation = 100)),
    "'grid' must contain a latitude column"
  )

  # a required variable is absent from the cell
  short_cols <- inputs$data
  short_cols[[1]]$temp_max <- NULL
  expect_error(.cp_call(inputs, data = short_cols), "missing columns: temp_max")

  # the cell has the wrong number of rows for the date vector
  short_rows <- inputs$data
  short_rows[[1]] <- short_rows[[1]][-1, ]
  expect_error(.cp_call(inputs, data = short_rows), "row count does not match")
})

test_that("apply_climate_perturbations: a 'y' column is accepted in place of 'lat'", {
  inputs <- .cp_valid_inputs()

  # the latitude column may be named either way; this must not error
  expect_no_error(
    .cp_call(
      inputs,
      grid = data.frame(id = 1, y = 0),
      temp_transient = FALSE,
      precip_transient = FALSE,
      precip_occurrence_transient = FALSE,
      compute_pet = FALSE,
      scale_var_with_mean = FALSE,
      enforce_target_mean = FALSE,
      diagnostic = FALSE,
      verbose = FALSE
    )
  )
})

test_that("apply_climate_perturbations: safety-rail arguments are validated", {
  inputs <- .cp_valid_inputs()

  expect_error(
    .cp_call(inputs, precip_intensity_threshold = -1),
    "'precip_intensity_threshold' must be a single finite numeric >= 0"
  )
  expect_error(
    .cp_call(inputs, precip_intensity_threshold = c(0, 1)),
    "'precip_intensity_threshold' must be a single finite numeric >= 0"
  )
  expect_error(
    .cp_call(inputs, precip_cap_mm_day = 0),
    "'precip_cap_mm_day' must be a single positive finite numeric"
  )
  expect_error(
    .cp_call(inputs, precip_cap_mm_day = Inf),
    "'precip_cap_mm_day' must be a single positive finite numeric"
  )
  expect_error(
    .cp_call(inputs, precip_floor_mm_day = -0.5),
    "'precip_floor_mm_day' must be a single finite numeric >= 0"
  )
  expect_error(
    .cp_call(inputs, precip_cap_quantile = 1),
    "'precip_cap_quantile' must be a single numeric in \\(0, 1\\)"
  )
  expect_error(
    .cp_call(inputs, precip_cap_quantile = 0),
    "'precip_cap_quantile' must be a single numeric in \\(0, 1\\)"
  )
  expect_error(
    .cp_call(inputs, seed = "abc"),
    "'seed' must be a single finite numeric/integer"
  )
})

test_that("apply_climate_perturbations: matrix change factors are validated against the date span", {
  n_years <- 2L
  inputs <- .cp_valid_inputs(n_years)

  # wrong number of month columns
  expect_error(
    .cp_call(inputs, precip_mean_factor = matrix(1, nrow = n_years, ncol = 11L)),
    "'precip_mean_factor' must have 12 columns \\(months\\)"
  )

  # right shape, wrong number of years
  expect_error(
    .cp_call(inputs, precip_mean_factor = matrix(1, nrow = n_years + 1L, ncol = 12L)),
    "'precip_mean_factor' must have nrow == number of years in 'date'"
  )

  # correct shape but containing a non-finite entry
  bad <- matrix(1, nrow = n_years, ncol = 12L)
  bad[1, 1] <- NA_real_
  expect_error(
    .cp_call(inputs, precip_mean_factor = bad),
    "'precip_mean_factor' contains non-finite values"
  )

  # a well-formed matrix is accepted where a length-12 vector would be
  expect_no_error(
    .cp_call(
      inputs,
      precip_mean_factor = matrix(1, nrow = n_years, ncol = 12L),
      temp_transient = FALSE,
      precip_transient = FALSE,
      precip_occurrence_transient = FALSE,
      compute_pet = FALSE,
      scale_var_with_mean = FALSE,
      enforce_target_mean = FALSE,
      diagnostic = FALSE,
      verbose = FALSE
    )
  )
})


testthat::test_that("apply_climate_perturbations restores the caller's RNG state", {
  grid <- data.frame(lat = c(10, 11))
  dates <- generate_noleap_dates(as.Date("2020-01-01"), 365)
  cell <- data.frame(
    precip = pmax(stats::rnorm(365, 3, 3), 0),
    temp = stats::rnorm(365, 20, 5),
    temp_min = stats::rnorm(365, 15, 5),
    temp_max = stats::rnorm(365, 25, 5)
  )
  dat <- list(cell, cell)

  run <- function() {
    apply_climate_perturbations(
      data = dat, grid = grid, date = dates,
      precip_mean_factor = rep(0.9, 12),
      temp_delta = rep(1, 12),
      seed = 42L, verbose = FALSE, diagnostic = FALSE,
      compute_pet = FALSE
    )
  }

  # Seeding inside the function must not move the caller's stream.
  set.seed(99)
  before <- .Random.seed
  invisible(run())
  testthat::expect_identical(.Random.seed, before)

  # And the caller's next draw is unaffected.
  set.seed(99)
  expected <- stats::runif(3)
  set.seed(99)
  invisible(run())
  testthat::expect_equal(stats::runif(3), expected)
})

testthat::test_that("apply_climate_perturbations names its own arguments in errors", {
  grid <- data.frame(lat = 10)
  dates <- generate_noleap_dates(as.Date("2020-01-01"), 365)

  testthat::expect_error(
    apply_climate_perturbations(data = NULL, grid = grid, date = dates,
                                precip_mean_factor = rep(1, 12),
                                temp_delta = rep(0, 12)),
    "'data' must not be NULL", fixed = TRUE
  )
  testthat::expect_error(
    apply_climate_perturbations(data = list(), grid = grid, date = NULL,
                                precip_mean_factor = rep(1, 12),
                                temp_delta = rep(0, 12)),
    "'date' must not be NULL", fixed = TRUE
  )

  # The latitude message names both accepted column names.
  bad_grid <- data.frame(z = 10)
  testthat::expect_error(
    apply_climate_perturbations(data = list(data.frame()), grid = bad_grid,
                                date = dates, precip_mean_factor = rep(1, 12),
                                temp_delta = rep(0, 12)),
    "'lat' or 'y'", fixed = TRUE
  )
})


testthat::test_that("temp_range_factor scales the diurnal range about its midpoint", {
  grid <- data.frame(lat = c(10, 11))
  dates <- generate_noleap_dates(as.Date("2020-01-01"), 730)
  set.seed(8)
  n <- length(dates)
  cell <- data.frame(
    precip = pmax(stats::rnorm(n, 3, 3), 0),
    temp = stats::rnorm(n, 25, 3),
    temp_min = stats::rnorm(n, 20, 3),
    temp_max = stats::rnorm(n, 30, 3)
  )
  cell$temp_max <- pmax(cell$temp_max, cell$temp_min + 1)
  dat <- list(cell, cell)

  run <- function(...) {
    apply_climate_perturbations(
      data = dat, grid = grid, date = dates,
      precip_mean_factor = rep(1, 12), temp_delta = rep(0, 12),
      temp_transient = FALSE, verbose = FALSE, diagnostic = FALSE,
      compute_pet = TRUE, ...
    )
  }

  base <- run()[[1]]
  wide <- run(temp_range_factor = rep(1.2, 12))[[1]]

  base_dtr <- base$temp_max - base$temp_min
  wide_dtr <- wide$temp_max - wide$temp_min

  testthat::expect_equal(wide_dtr, base_dtr * 1.2)

  # The midpoint of min and max is preserved, so the daily mean is untouched.
  testthat::expect_equal((wide$temp_max + wide$temp_min) / 2,
                         (base$temp_max + base$temp_min) / 2)
  testthat::expect_equal(wide$temp, base$temp)

  # PET goes as sqrt(range), so a 1.2x range is a sqrt(1.2)x PET.
  testthat::expect_equal(mean(wide$pet) / mean(base$pet), sqrt(1.2), tolerance = 1e-6)

  # NULL is the previous behaviour exactly.
  testthat::expect_equal(run(temp_range_factor = NULL)[[1]], base)

  # A shrinking range lowers PET, and the range stays positive.
  narrow <- run(temp_range_factor = rep(0.5, 12))[[1]]
  testthat::expect_true(all(narrow$temp_max >= narrow$temp_min))
  testthat::expect_lt(mean(narrow$pet), mean(base$pet))
})

# ==============================================================================
# Stage handoff: prepare_evaluation_data() -> apply_climate_perturbations()
# ==============================================================================

# The README and the Climate Perturbations vignette document perturbation as a
# stage applied to an already-generated series, chained through
# prepare_evaluation_data(). That chain is a documented workflow with no
# wrapper function holding it together, so these tests pin the interface
# instead. generate_weather() is not called -- a stub gen_output carries the
# same contract at a fraction of the cost.

.make_gen_output_stub <- function(obs_dates, n_sim_days, n_realizations = 2L,
                                  seed = 11L) {
  set.seed(seed)
  sim_dates <- .make_dates_noleap_years("2100-01-01", n_sim_days / 365L)

  resampled <- as.data.frame(
    lapply(
      seq_len(n_realizations),
      function(i) obs_dates[sample.int(length(obs_dates), length(sim_dates),
                                       replace = TRUE)]
    )
  )
  names(resampled) <- paste0("rlz_", seq_len(n_realizations))

  list(resampled = resampled, dates = sim_dates)
}

.make_obs_for_handoff <- function(date, n_cells = 2L) {
  lapply(seq_len(n_cells), function(i) {
    n <- length(date)
    set.seed(100L + i)
    cell <- data.frame(
      precip   = .make_precip_with_min_wet(date, min_wet = 15L, seed = 100L + i),
      temp     = stats::rnorm(n, 25, 3),
      temp_min = stats::rnorm(n, 20, 3)
    )
    cell$temp_max <- pmax(cell$temp_min + 5, stats::rnorm(n, 30, 3))
    cell
  })
}

test_that("prepare_evaluation_data() output feeds apply_climate_perturbations() unchanged", {
  vars <- c("precip", "temp", "temp_min", "temp_max")
  obs_dates <- .make_dates_noleap_years("2000-01-01", 4L)
  obs_data <- .make_obs_for_handoff(obs_dates, n_cells = 2L)
  grid <- data.frame(id = 1:2, lat = c(0, 10))

  gen_output <- .make_gen_output_stub(obs_dates, n_sim_days = 365L * 3L)

  eval_data <- prepare_evaluation_data(
    gen_output = gen_output,
    obs_data   = obs_data,
    obs_dates  = obs_dates,
    grid_ids   = seq_along(obs_data),
    variables  = vars,
    verbose    = FALSE
  )

  sim <- eval_data$sim_data[[1]]

  # The documented shape: one data.frame per grid cell, date column first.
  expect_length(sim, length(obs_data))
  expect_equal(names(sim[[1]]), c("date", vars))

  # The extra 'date' column must not trip the required-columns check, which is
  # a setdiff() and so tolerates columns it did not ask for.
  out <- apply_climate_perturbations(
    data               = sim,
    grid               = grid,
    date               = sim[[1]]$date,
    precip_mean_factor = rep(0.7, 12),
    precip_var_factor  = rep(1.0, 12),
    temp_delta         = rep(2, 12),
    temp_transient     = FALSE,
    precip_transient   = FALSE,
    compute_pet        = FALSE,
    seed               = 42L,
    diagnostic         = FALSE,
    verbose            = FALSE
  )

  expect_length(out, length(sim))
  expect_equal(names(out[[1]]), c("date", vars))
  expect_equal(nrow(out[[1]]), nrow(sim[[1]]))

  # The perturbation lands on target through the chain.
  expect_equal(mean(out[[1]]$temp) - mean(sim[[1]]$temp), 2, tolerance = 1e-8)
  expect_equal(mean(out[[1]]$precip) / mean(sim[[1]]$precip), 0.7,
               tolerance = 0.05)
})

test_that("the date vector must come from sim_data, not gen_output$dates", {
  # prepare_evaluation_data() drops incomplete years, so its row count can be
  # shorter than gen_output$dates. Taking the date vector from the returned
  # frames is what keeps the two aligned; this pins the failure a naive recipe
  # would hit rather than leaving it for a user to discover.
  vars <- c("precip", "temp", "temp_min", "temp_max")
  obs_dates <- .make_dates_noleap_years("2000-01-01", 4L)
  obs_data <- .make_obs_for_handoff(obs_dates, n_cells = 2L)
  grid <- data.frame(id = 1:2, lat = c(0, 10))

  gen_output <- .make_gen_output_stub(obs_dates, n_sim_days = 365L * 3L)

  # A water-year boundary drops the partial first and last years.
  eval_data <- prepare_evaluation_data(
    gen_output       = gen_output,
    obs_data         = obs_data,
    obs_dates        = obs_dates,
    grid_ids         = seq_along(obs_data),
    variables        = vars,
    year_start_month = 10L,
    verbose          = FALSE
  )

  sim <- eval_data$sim_data[[1]]
  expect_lt(nrow(sim[[1]]), length(gen_output$dates))

  # The trimmed vector works.
  expect_no_error(
    apply_climate_perturbations(
      data = sim, grid = grid, date = sim[[1]]$date,
      precip_mean_factor = rep(0.7, 12), temp_delta = rep(2, 12),
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE
    )
  )

  # The untrimmed one does not, and says so by row count.
  expect_error(
    apply_climate_perturbations(
      data = sim, grid = grid, date = gen_output$dates,
      precip_mean_factor = rep(0.7, 12), temp_delta = rep(2, 12),
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE
    ),
    "row count"
  )
})

test_that("matrix factors are sized by calendar years spanned, not generator n_years", {
  # year_idx is calendar-based (cal_year - min(cal_year) + 1), so a series
  # covering N water years spans N + 1 calendar years and needs (N + 1) x 12.
  vars <- c("precip", "temp", "temp_min", "temp_max")
  obs_dates <- .make_dates_noleap_years("2000-01-01", 4L)
  obs_data <- .make_obs_for_handoff(obs_dates, n_cells = 2L)
  grid <- data.frame(id = 1:2, lat = c(0, 10))

  gen_output <- .make_gen_output_stub(obs_dates, n_sim_days = 365L * 3L)

  eval_data <- prepare_evaluation_data(
    gen_output = gen_output, obs_data = obs_data, obs_dates = obs_dates,
    grid_ids = seq_along(obs_data), variables = vars,
    year_start_month = 10L, verbose = FALSE
  )
  sim <- eval_data$sim_data[[1]]
  sim_date <- sim[[1]]$date

  n_cal <- length(unique(format(sim_date, "%Y")))
  n_wy  <- length(unique(compute_water_year(sim_date, 10L)))
  expect_gt(n_cal, n_wy)

  run_with_nrow <- function(nr) {
    apply_climate_perturbations(
      data = sim, grid = grid, date = sim_date,
      precip_mean_factor = matrix(0.7, nrow = nr, ncol = 12),
      precip_var_factor  = matrix(1.0, nrow = nr, ncol = 12),
      temp_delta = rep(2, 12), compute_pet = FALSE,
      diagnostic = FALSE, verbose = FALSE
    )
  }

  expect_no_error(run_with_nrow(n_cal))
  expect_error(run_with_nrow(n_wy), "nrow")

  # The length-12 vector form is indifferent to the distinction.
  expect_no_error(
    apply_climate_perturbations(
      data = sim, grid = grid, date = sim_date,
      precip_mean_factor = rep(0.7, 12), temp_delta = rep(2, 12),
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE
    )
  )
})

# ==============================================================================
# Transient factors reject year-varying matrices
# ==============================================================================

# A transient factor is specified by its end state and ramps from 1 to 2f-1,
# reading row 1 alone. Combining that with a year-varying matrix used to discard
# rows 2:n silently -- a matrix encoding a 50 percent reduction came back as no
# change at all.

test_that("a year-varying matrix with transient = TRUE is an error, not a silent discard", {
  n_years <- 3L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  dat <- list(.handoff_cell(date))
  grid <- data.frame(id = 1L, lat = 0)

  varying <- matrix(0.5, nrow = n_years, ncol = 12)
  varying[1, ] <- 1.0

  run <- function(...) {
    apply_climate_perturbations(
      data = dat, grid = grid, date = date,
      temp_delta = rep(0, 12), compute_pet = FALSE,
      diagnostic = FALSE, verbose = FALSE, seed = 1L, ...
    )
  }

  expect_error(
    run(precip_mean_factor = varying, precip_transient = TRUE),
    "precip_mean_factor.*varies by year"
  )

  # Same matrix is fine when the flag says apply it year by year.
  expect_no_error(run(precip_mean_factor = varying, precip_transient = FALSE))

  # A length-12 vector is the transient specification and stays accepted.
  expect_no_error(run(precip_mean_factor = rep(0.5, 12), precip_transient = TRUE))

  # A constant-row matrix is equivalent to a vector, so it is still accepted.
  expect_no_error(
    run(precip_mean_factor = matrix(0.5, nrow = n_years, ncol = 12),
        precip_transient = TRUE)
  )
})

test_that("the transient guard names the argument the caller actually supplied", {
  n_years <- 3L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  dat <- list(.handoff_cell(date))
  grid <- data.frame(id = 1L, lat = 0)

  varying <- matrix(0.5, nrow = n_years, ncol = 12)
  varying[1, ] <- 1.0
  flat <- rep(1, 12)

  # scale_var_with_mean = TRUE derives the variance factor from the mean, so a
  # varying mean must be reported against precip_mean_factor rather than against
  # a precip_var_factor the caller never passed.
  expect_error(
    apply_climate_perturbations(
      data = dat, grid = grid, date = date,
      precip_mean_factor = varying, scale_var_with_mean = TRUE,
      temp_delta = rep(0, 12), precip_transient = TRUE,
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE
    ),
    "precip_mean_factor"
  )

  # With scale_var_with_mean = FALSE the variance factor is the caller's own.
  expect_error(
    apply_climate_perturbations(
      data = dat, grid = grid, date = date,
      precip_mean_factor = flat, precip_var_factor = varying,
      scale_var_with_mean = FALSE,
      temp_delta = rep(0, 12), precip_transient = TRUE,
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE
    ),
    "precip_var_factor.*varies by year"
  )

  # Occurrence and diurnal range route through the same guard.
  expect_error(
    apply_climate_perturbations(
      data = dat, grid = grid, date = date,
      precip_mean_factor = flat, precip_occurrence_factor = varying,
      temp_delta = rep(0, 12), precip_occurrence_transient = TRUE,
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE
    ),
    "precip_occurrence_factor.*varies by year"
  )

  expect_error(
    apply_climate_perturbations(
      data = dat, grid = grid, date = date,
      precip_mean_factor = flat, temp_range_factor = varying,
      temp_delta = rep(0, 12), temp_transient = TRUE,
      compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE
    ),
    "temp_range_factor.*varies by year"
  )
})

test_that("a transient ramp averages to the requested end-state factor", {
  # The documented property: ramp from 1 to 2f-1, so the period mean is f.
  n_years <- 5L
  date <- .make_dates_noleap_years("2000-01-01", n_years)
  dat <- list(.handoff_cell(date))
  grid <- data.frame(id = 1L, lat = 0)

  out <- apply_climate_perturbations(
    data = dat, grid = grid, date = date,
    precip_mean_factor = rep(0.6, 12),
    temp_delta = rep(0, 12), precip_transient = TRUE,
    compute_pet = FALSE, diagnostic = FALSE, verbose = FALSE, seed = 1L
  )

  yidx <- as.integer(format(date, "%Y")) - 2000L + 1L
  obs_ann <- tapply(dat[[1]]$precip, yidx, sum)
  new_ann <- tapply(out[[1]]$precip, yidx, sum)
  ratio <- as.numeric(new_ann / obs_ann)

  # Declines across the period, and lands near 0.6 on average.
  expect_lt(ratio[n_years], ratio[1])
  expect_equal(mean(ratio), 0.6, tolerance = 0.15)
})
