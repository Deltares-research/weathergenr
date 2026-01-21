# Functions tested (relative paths):
# - R/generator.R: generate_weather()

make_generator_inputs <- function() {
  set.seed(123)
  obs_dates <- seq.Date(as.Date("2000-01-01"), as.Date("2015-12-31"), by = "day")

  n <- length(obs_dates)
  precip <- 5 + sin(seq(0, 20 * pi, length.out = n)) + rnorm(n, sd = 0.5)
  temp <- 10 + cos(seq(0, 10 * pi, length.out = n)) + rnorm(n, sd = 0.3)

  obs_data <- list(data.frame(precip = precip, temp = temp))
  obs_grid <- data.frame(xind = 1, yind = 1, x = 0, y = 0)

  list(
    obs_data = obs_data,
    obs_grid = obs_grid,
    obs_dates = obs_dates,
    vars = c("precip", "temp")
  )
}

testthat::test_that("generate_weather returns expected structure quickly", {
  inputs <- make_generator_inputs()
  out_dir <- file.path(tempdir(), "weathergenr-generator")

  out <- suppressWarnings(generate_weather(
    obs_data = inputs$obs_data,
    obs_grid = inputs$obs_grid,
    obs_dates = inputs$obs_dates,
    vars = inputs$vars,
    n_years = 16,
    start_year = 2020,
    year_start_month = 1,
    n_realizations = 1,
    warm_var = "precip",
    warm_signif = 0.8,
    warm_pool_size = 3,
    annual_knn_n = 3,
    wet_q = 0.3,
    extreme_q = 0.8,
    out_dir = out_dir,
    seed = 42,
    parallel = FALSE,
    verbose = FALSE
  ))

  testthat::expect_true(all(c("resampled", "dates") %in% names(out)))
  testthat::expect_s3_class(out$dates, "Date")
  testthat::expect_length(out$dates, 16 * 365)
  testthat::expect_true(is.data.frame(out$resampled))
  testthat::expect_equal(nrow(out$resampled), length(out$dates))
  testthat::expect_true(all(grepl("^rlz_", names(out$resampled))))
  testthat::expect_false(anyNA(out$resampled))
})

testthat::test_that("generate_weather validates wet/extreme thresholds", {
  inputs <- make_generator_inputs()

  testthat::expect_error(
    generate_weather(
      obs_data = inputs$obs_data,
      obs_grid = inputs$obs_grid,
      obs_dates = inputs$obs_dates,
      vars = inputs$vars,
      n_years = 16,
      start_year = 2020,
      year_start_month = 1,
      n_realizations = 1,
      warm_var = "precip",
      warm_signif = 0.8,
      warm_pool_size = 3,
      annual_knn_n = 3,
      wet_q = 0.8,
      extreme_q = 0.8,
      out_dir = tempdir(),
      seed = 42,
      parallel = FALSE,
      verbose = FALSE
    ),
    "extreme_q must be greater than wet_q"
  )
})
