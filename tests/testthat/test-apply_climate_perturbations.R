testthat::test_that("apply_climate_perturbations: returns list-of-data.frames with same structure", {
  n_years <- 2L
  n_days <- 365L * n_years
  date <- as.Date("2000-01-01") + 0:(n_days - 1L)

  grid <- data.frame(id = 1:2, lat = c(0, 10))

  # deterministic synthetic series (no NA)
  mk_cell <- function(seed) {
    set.seed(seed)
    month <- as.integer(format(date, "%m"))
    prcp <- rgamma(n_days, shape = 1.2, scale = 5)
    prcp[sample.int(n_days, 150)] <- 0
    temp <- 20 + sin(2 * pi * (1:n_days) / 365) * 3
    out <- data.frame(
      prcp = prcp,
      temp = temp,
      temp_min = temp - 2,
      temp_max = temp + 2
    )
    out
  }

  data <- list(mk_cell(1), mk_cell(2))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    prcp_mean_factor = rep(1, 12),
    prcp_var_factor = rep(1, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    prcp_transient = FALSE,
    prcp_occurrence_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,   # important: avoid rescaling in QM
    seed = 1,
    verbose = FALSE
  )

  testthat::expect_true(is.list(out))
  testthat::expect_length(out, length(data))

  for (i in seq_along(out)) {
    testthat::expect_true(is.data.frame(out[[i]]))
    testthat::expect_equal(nrow(out[[i]]), n_days)
    testthat::expect_true(all(c("prcp", "temp", "temp_min", "temp_max") %in% names(out[[i]])))
  }
})

testthat::test_that("apply_climate_perturbations: temperature deltas apply (step change)", {
  n_years <- 2L
  n_days <- 365L * n_years
  date <- as.Date("2000-01-01") + 0:(n_days - 1L)

  grid <- data.frame(id = 1, lat = 0)

  prcp <- rep(1, n_days)
  temp <- rep(10, n_days)

  data <- list(data.frame(
    prcp = prcp,
    temp = temp,
    temp_min = temp - 1,
    temp_max = temp + 1
  ))

  delta <- rep(2, 12)

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    prcp_mean_factor = rep(1, 12),
    prcp_var_factor = rep(1, 12),
    temp_delta = delta,
    temp_transient = FALSE,
    prcp_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    verbose = FALSE
  )[[1]]

  testthat::expect_equal(out$temp, temp + 2)
  testthat::expect_equal(out$temp_min, (temp - 1) + 2)
  testthat::expect_equal(out$temp_max, (temp + 1) + 2)
})

testthat::test_that("apply_climate_perturbations: transient temperature ramp has 0 at first year and ~2*delta at last year", {
  n_years <- 3L
  n_days <- 365L * n_years
  date <- as.Date("2000-01-01") + 0:(n_days - 1L)

  grid <- data.frame(id = 1, lat = 0)

  temp <- rep(0, n_days)
  data <- list(data.frame(
    prcp = rep(1, n_days),
    temp = temp,
    temp_min = temp,
    temp_max = temp
  ))

  delta <- rep(1, 12)

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    prcp_mean_factor = rep(1, 12),
    prcp_var_factor = rep(1, 12),
    temp_delta = delta,
    temp_transient = TRUE,
    prcp_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    verbose = FALSE
  )[[1]]

  cal_year <- as.integer(format(date, "%Y"))
  year_idx <- cal_year - min(cal_year) + 1L

  # Pick a day in January of year 1 and year 3
  i1 <- which(year_idx == 1 & format(date, "%m") == "01")[1]
  i3 <- which(year_idx == 3 & format(date, "%m") == "01")[1]

  testthat::expect_equal(out$temp[i1], 0, tolerance = 1e-12)
  testthat::expect_equal(out$temp[i3], 2, tolerance = 1e-12)  # 2*delta, delta=1
})

testthat::test_that("apply_climate_perturbations: precipitation unchanged under identity factors when enforce_target_mean=FALSE", {
  n_years <- 2L
  n_days <- 365L * n_years
  date <- as.Date("2000-01-01") + 0:(n_days - 1L)

  grid <- data.frame(id = 1, lat = 0)

  set.seed(3)
  prcp <- rgamma(n_days, shape = 1.5, scale = 4)
  prcp[sample.int(n_days, 120)] <- 0

  data <- list(data.frame(
    prcp = prcp,
    temp = rep(0, n_days),
    temp_min = rep(0, n_days),
    temp_max = rep(0, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    prcp_mean_factor = rep(1, 12),
    prcp_var_factor = rep(1, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    prcp_transient = FALSE,
    prcp_occurrence_factor = NULL,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,  # key for identity-like behavior
    seed = 1,
    verbose = FALSE
  )[[1]]

  wet <- prcp > 0
  testthat::expect_identical(out$prcp[!wet], prcp[!wet])
  testthat::expect_lt(mean(abs(out$prcp[wet] - prcp[wet])), 1e-6)
})

testthat::test_that("apply_climate_perturbations: occurrence factor changes wet-day frequency in expected direction", {
  n_years <- 2L
  n_days <- 365L * n_years
  date <- as.Date("2000-01-01") + 0:(n_days - 1L)
  month <- as.integer(format(date, "%m"))

  grid <- data.frame(id = 1, lat = 0)

  set.seed(7)
  prcp <- rgamma(n_days, shape = 1.2, scale = 4)
  # enforce many dry days
  prcp[sample.int(n_days, 250)] <- 0

  data <- list(data.frame(
    prcp = prcp,
    temp = rep(0, n_days),
    temp_min = rep(0, n_days),
    temp_max = rep(0, n_days)
  ))

  # no intensity change, only occurrence
  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    prcp_mean_factor = rep(1, 12),
    prcp_var_factor = rep(1, 12),
    prcp_occurrence_factor = rep(1.5, 12),  # increase wet frequency
    prcp_intensity_threshold = 0,
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    prcp_transient = FALSE,
    prcp_occurrence_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    seed = 42,
    verbose = FALSE
  )[[1]]

  wet0 <- prcp > 0
  wet1 <- out$prcp > 0

  # Expect wet-day count increases overall (stochastic rounding; allow ties rarely)
  testthat::expect_gte(sum(wet1), sum(wet0))

  # Added wet days should have positive amounts
  added <- (!wet0) & wet1
  if (any(added)) {
    testthat::expect_true(all(out$prcp[added] > 0))
  }
})

testthat::test_that("apply_climate_perturbations: safety rails enforce cap and floor", {
  n_years <- 2L
  n_days <- 365L * n_years
  date <- as.Date("2000-01-01") + 0:(n_days - 1L)

  grid <- data.frame(id = 1, lat = 0)

  set.seed(9)
  prcp <- rgamma(n_days, shape = 0.8, scale = 20)
  prcp[sample.int(n_days, 100)] <- 0

  data <- list(data.frame(
    prcp = prcp,
    temp = rep(0, n_days),
    temp_min = rep(0, n_days),
    temp_max = rep(0, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    prcp_mean_factor = rep(1.2, 12),
    prcp_var_factor = rep(1.0, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    prcp_transient = FALSE,
    compute_pet = FALSE,
    scale_var_with_mean = FALSE,
    exaggerate_extremes = TRUE,
    extreme_k = 1.5,
    enforce_target_mean = TRUE,
    prcp_floor_mm_day = 0.2,
    prcp_cap_mm_day = 50,
    seed = 1,
    verbose = FALSE
  )[[1]]

  wet <- out$prcp > 0
  if (any(wet)) {
    testthat::expect_true(all(out$prcp[wet] >= 0.2))
    testthat::expect_true(max(out$prcp, na.rm = TRUE) <= 50)
  }
})

testthat::test_that("apply_climate_perturbations: PET is added when compute_pet=TRUE", {
  testthat::skip_if_not(exists("calculate_monthly_pet", where = asNamespace("weathergenr"), inherits = FALSE))

  n_years <- 2L
  n_days <- 365L * n_years
  date <- as.Date("2000-01-01") + 0:(n_days - 1L)

  grid <- data.frame(id = 1, lat = 10)

  data <- list(data.frame(
    prcp = rep(1, n_days),
    temp = rep(20, n_days),
    temp_min = rep(18, n_days),
    temp_max = rep(22, n_days)
  ))

  out <- apply_climate_perturbations(
    data = data,
    grid = grid,
    date = date,
    prcp_mean_factor = rep(1, 12),
    prcp_var_factor = rep(1, 12),
    temp_delta = rep(0, 12),
    temp_transient = FALSE,
    prcp_transient = FALSE,
    compute_pet = TRUE,
    pet_method = "hargreaves",
    scale_var_with_mean = FALSE,
    exaggerate_extremes = FALSE,
    enforce_target_mean = FALSE,
    verbose = FALSE
  )[[1]]

  testthat::expect_true("pet" %in% names(out))
  testthat::expect_length(out$pet, n_days)
  testthat::expect_true(all(is.finite(out$pet)))
  testthat::expect_true(all(out$pet >= 0))
})
