# Functions tested (relative paths):
# - R/generator.R: generate_weather()
# - R/generator.R: run_weather_generator()

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

testthat::test_that("generate_weather validates required daily variables", {
  inputs <- make_generator_inputs()

  testthat::expect_error(
    generate_weather(
      obs_data = inputs$obs_data,
      obs_grid = inputs$obs_grid,
      obs_dates = inputs$obs_dates,
      vars = "precip",
      n_years = 2,
      start_year = 2020,
      year_start_month = 1,
      n_realizations = 1,
      warm_var = "precip",
      warm_signif = 0.8,
      warm_pool_size = 3,
      annual_knn_n = 3,
      wet_q = 0.3,
      extreme_q = 0.8,
      out_dir = tempdir(),
      seed = 42,
      parallel = FALSE,
      verbose = FALSE
    ),
    "vars must include 'precip' and 'temp' for daily disaggregation"
  )
})

# ==============================================================================
# run_weather_generator(): top-level orchestration
#
# The three downstream steps are mocked throughout. The point of these tests is
# the wiring -- validation, grid-ID derivation, argument threading, logging and
# sink hygiene -- not the numerics of generation, which the generate_weather()
# tests above and test-warm.R / test-resample.R already cover. Mocking also
# keeps this whole block well under a second.
# ==============================================================================

# Minimal well-formed stubs for the three steps run_weather_generator() calls.
# Each records the arguments it received into `calls` so a test can assert on
# what was threaded through.
.mock_run_steps <- function(calls) {
  sim_dates <- seq.Date(as.Date("2021-01-01"), by = "day", length.out = 4L)

  list(
    generate_weather = function(...) {
      calls$generate <- list(...)
      list(
        resampled = data.frame(rlz_1 = sim_dates),
        dates = sim_dates
      )
    },
    prepare_evaluation_data = function(...) {
      calls$prepare <- list(...)
      list(sim_data = list(), obs_data = list())
    },
    evaluate_weather_generator = function(...) {
      calls$evaluate <- list(...)
      structure(list(mocked = TRUE), class = "weather_assessment")
    }
  )
}

# Applies the stubs for the duration of the calling test.
.local_mock_run_steps <- function(calls, frame = parent.frame()) {
  m <- .mock_run_steps(calls)
  testthat::local_mocked_bindings(
    generate_weather           = m$generate_weather,
    prepare_evaluation_data    = m$prepare_evaluation_data,
    evaluate_weather_generator = m$evaluate_weather_generator,
    .env                       = frame
  )
  invisible(m)
}

.run_config <- function(...) {
  utils::modifyList(
    list(
      vars = c("precip", "temp"),
      n_realizations = 1L,
      n_years = 2,
      start_year = 2021,
      year_start_month = 1,
      warm_var = "precip",
      wet_q = 0.3,
      extreme_q = 0.8,
      parallel = FALSE,
      seed = 42,
      verbose = FALSE,
      save_plots = FALSE
    ),
    list(...)
  )
}

testthat::test_that("run_weather_generator rejects a malformed config", {
  inputs <- make_generator_inputs()

  args <- list(
    obs_data  = inputs$obs_data,
    obs_grid  = inputs$obs_grid,
    obs_dates = inputs$obs_dates,
    out_dir   = tempdir()
  )

  testthat::expect_error(
    do.call(run_weather_generator, c(args, list(config = "not-a-list"))),
    "'config' must be a list"
  )

  testthat::expect_error(
    do.call(run_weather_generator, c(args, list(config = list(n_realizations = 1L)))),
    "'config\\$vars' must be provided"
  )

  testthat::expect_error(
    do.call(
      run_weather_generator,
      c(args, list(config = list(vars = character(0), n_realizations = 1L)))
    ),
    "'config\\$vars' must be provided"
  )

  testthat::expect_error(
    do.call(
      run_weather_generator,
      c(args, list(config = list(vars = c("precip", "temp"))))
    ),
    "'config\\$n_realizations' must be provided"
  )
})

testthat::test_that("run_weather_generator returns generation and evaluation output", {
  inputs <- make_generator_inputs()
  calls <- new.env(parent = emptyenv())
  .local_mock_run_steps(calls)

  out <- run_weather_generator(
    obs_data  = inputs$obs_data,
    obs_grid  = inputs$obs_grid,
    obs_dates = inputs$obs_dates,
    out_dir   = tempdir(),
    config    = .run_config()
  )

  testthat::expect_true(all(c("resampled", "dates") %in% names(out$gen_output)))
  testthat::expect_s3_class(out$evaluation, "weather_assessment")

  # log_messages defaults to FALSE, so no log is written and no sink is left open.
  #
  # Note the return shape here: `res$log_path <- NULL` *removes* the element
  # rather than storing a NULL, so with logging off the returned list carries
  # only two components. `out$log_path` is still NULL -- which is the contract
  # callers actually rely on -- but `names(out)` does not list it, unlike the
  # @return block, which describes log_path as a component that is NULL. This
  # test pins the behavior as it stands rather than asserting the docs.
  testthat::expect_named(out, c("gen_output", "evaluation"))
  testthat::expect_null(out$log_path)
  testthat::expect_identical(sink.number(), 0L)
})

testthat::test_that("run_weather_generator threads config and out_dir into generate_weather", {
  inputs <- make_generator_inputs()
  calls <- new.env(parent = emptyenv())
  .local_mock_run_steps(calls)

  out_dir <- file.path(tempdir(), "weathergenr-run-threading")
  config <- .run_config(n_realizations = 3L, seed = 7, wet_q = 0.25, extreme_q = 0.95)

  run_weather_generator(
    obs_data  = inputs$obs_data,
    obs_grid  = inputs$obs_grid,
    obs_dates = inputs$obs_dates,
    out_dir   = out_dir,
    config    = config
  )

  gen <- calls$generate
  testthat::expect_identical(gen$vars, config$vars)
  testthat::expect_identical(gen$n_realizations, 3L)
  testthat::expect_identical(gen$seed, 7)
  testthat::expect_identical(gen$wet_q, 0.25)
  testthat::expect_identical(gen$extreme_q, 0.95)
  testthat::expect_identical(gen$warm_var, "precip")
  testthat::expect_false(gen$parallel)

  # out_dir is the function argument, not a config field
  testthat::expect_identical(gen$out_dir, out_dir)

  # eval_max_grids is likewise a function argument and reaches the evaluator
  testthat::expect_identical(calls$evaluate$eval_max_grids, 25L)
  testthat::expect_identical(calls$evaluate$n_realizations, 3L)
})

testthat::test_that("run_weather_generator derives grid IDs from each supported obs_grid shape", {
  inputs <- make_generator_inputs()

  grid_id_for <- function(obs_grid) {
    calls <- new.env(parent = emptyenv())
    .local_mock_run_steps(calls)
    run_weather_generator(
      obs_data  = inputs$obs_data,
      obs_grid  = obs_grid,
      obs_dates = inputs$obs_dates,
      out_dir   = tempdir(),
      config    = .run_config()
    )
    calls$prepare$grid_ids
  }

  # explicit ids win, from a list or a data.frame
  testthat::expect_identical(grid_id_for(list(id = c(11L, 12L))), c(11L, 12L))
  testthat::expect_identical(
    grid_id_for(data.frame(id = c(5L, 6L), x = 0, y = 0)),
    c(5L, 6L)
  )

  # otherwise fall back to a row/coordinate sequence
  testthat::expect_identical(grid_id_for(data.frame(x = c(0, 1, 2), y = 0)), 1:3)
  testthat::expect_identical(grid_id_for(list(lon = c(0, 1), lat = c(0, 1))), 1:2)
})

testthat::test_that("run_weather_generator errors when grid IDs cannot be derived", {
  inputs <- make_generator_inputs()
  calls <- new.env(parent = emptyenv())
  .local_mock_run_steps(calls)

  testthat::expect_error(
    run_weather_generator(
      obs_data  = inputs$obs_data,
      obs_grid  = list(crs = "EPSG:4326"),
      obs_dates = inputs$obs_dates,
      out_dir   = tempdir(),
      config    = .run_config()
    ),
    "Unable to derive grid IDs"
  )
})

testthat::test_that("run_weather_generator writes and merges a log when log_messages = TRUE", {
  inputs <- make_generator_inputs()
  calls <- new.env(parent = emptyenv())
  m <- .mock_run_steps(calls)

  # Emit a message and a warning from inside the run so the merge path is covered
  testthat::local_mocked_bindings(
    generate_weather = function(...) {
      message("mock-generator-message")
      warning("mock-generator-warning", call. = FALSE)
      m$generate_weather(...)
    },
    prepare_evaluation_data    = m$prepare_evaluation_data,
    evaluate_weather_generator = m$evaluate_weather_generator
  )

  # out_dir deliberately does not exist yet: run_weather_generator must create it
  out_dir <- file.path(tempdir(), "weathergenr-run-log")
  unlink(out_dir, recursive = TRUE)
  testthat::expect_false(dir.exists(out_dir))

  on.exit(
    {
      while (sink.number() > 0) sink(NULL)
      unlink(out_dir, recursive = TRUE)
    },
    add = TRUE
  )

  # Both are suppressed only at the console: run_weather_generator's calling
  # handlers deliberately do not muffle, so the message and warning still reach
  # the capture file before these outer suppressors swallow them.
  out <- suppressMessages(suppressWarnings(
    run_weather_generator(
      obs_data     = inputs$obs_data,
      obs_grid     = inputs$obs_grid,
      obs_dates    = inputs$obs_dates,
      out_dir      = out_dir,
      config       = .run_config(),
      log_messages = TRUE
    )
  ))

  testthat::expect_true(dir.exists(out_dir))
  testthat::expect_type(out$log_path, "character")
  testthat::expect_true(file.exists(out$log_path))
  testthat::expect_match(basename(out$log_path), "^log_[0-9]{8}_[0-9]{6}\\.txt$")

  log_lines <- readLines(out$log_path, warn = FALSE)
  testthat::expect_true(any(grepl("weathergenr run log", log_lines, fixed = TRUE)))

  # messages and warnings are merged in from the separate capture file...
  testthat::expect_true(any(grepl("MESSAGES/WARNINGS", log_lines, fixed = TRUE)))
  testthat::expect_true(any(grepl("mock-generator-message", log_lines, fixed = TRUE)))
  testthat::expect_true(any(grepl("mock-generator-warning", log_lines, fixed = TRUE)))

  # ...and the temporary capture file is cleaned up afterwards
  testthat::expect_length(list.files(out_dir, pattern = "_messages\\.txt$"), 0L)

  # the sink stack must be unwound even though logging was enabled
  testthat::expect_identical(sink.number(), 0L)
})
