# Functions tested (relative paths):
# - R/utils_internal.R: compute_skewness(), compute_kurtosis(), compute_spell_lengths(),
#   assess_moment_changes(), compute_area_averages(), .is_int_scalar(), .safe_cor(),
#   .format_pct(), .format_num(), .log(), format_elapsed()

testthat::test_that("compute_skewness handles NA, short input, and constant series", {
  testthat::expect_true(is.na(weathergenr:::compute_skewness(c(1, NA, 2))))
  testthat::expect_true(is.na(weathergenr:::compute_skewness(rep(5, 5))))

  x <- c(-1, 0, 1)
  testthat::expect_equal(weathergenr:::compute_skewness(x), 0)
})

testthat::test_that("compute_kurtosis handles NA, short input, and constant series", {
  testthat::expect_true(is.na(weathergenr:::compute_kurtosis(c(1, NA, 2, 3))))
  testthat::expect_true(is.na(weathergenr:::compute_kurtosis(rep(2, 6))))

  x <- c(-2, -1, 0, 1, 2)
  out <- weathergenr:::compute_kurtosis(x)
  testthat::expect_true(is.finite(out))
})

testthat::test_that("compute_spell_lengths returns expected run lengths", {
  x <- c(0, 0, 1, 2, 0, 0, 0, 3)
  dry <- weathergenr:::compute_spell_lengths(x, threshold = 0.5, below = TRUE)
  wet <- weathergenr:::compute_spell_lengths(x, threshold = 0.5, below = FALSE)

  testthat::expect_identical(dry, c(2L, 3L))
  testthat::expect_identical(wet, c(2L, 1L))

  none <- weathergenr:::compute_spell_lengths(c(0, 0, 0), threshold = 0.5, below = FALSE)
  testthat::expect_identical(none, numeric(0))
})

testthat::test_that("assess_moment_changes applies metric-specific thresholds", {
  moments_df <- data.frame(
    metric = c("mean", "variance", "cv", "sd", "skewness", "kurtosis"),
    pct_change = c(4, 10, 15, 14, 25, 31),
    stringsAsFactors = FALSE
  )

  out <- weathergenr:::assess_moment_changes(moments_df)
  testthat::expect_identical(out, c("excellent", "good", "good", "good", "acceptable", "poor"))
})

testthat::test_that("compute_area_averages returns expected structure", {
  grid1 <- data.frame(precip = c(1, 2, 3), temp = c(10, 11, 12))
  grid2 <- data.frame(precip = c(2, 3, 4), temp = c(9, 10, 11))
  obs_data <- list(g1 = grid1, g2 = grid2)
  wyear_idx <- 1:3
  wyear <- c(2001, 2001, 2002)
  vars <- c("precip", "temp")

  out <- compute_area_averages(obs_data, wyear_idx, wyear, vars)
  testthat::expect_true(all(c("daily", "annual") %in% names(out)))
  testthat::expect_true(all(c("precip", "temp", "wyear") %in% names(out$daily)))
  testthat::expect_true(all(c("precip", "temp", "wyear") %in% names(out$annual)))

  testthat::expect_equal(out$daily$precip, c(1.5, 2.5, 3.5))
  testthat::expect_equal(out$daily$temp, c(9.5, 10.5, 11.5))
  testthat::expect_equal(out$annual$precip, c(2, 3.5))
})

testthat::test_that(".is_int_scalar validates integer-valued scalars", {
  testthat::expect_true(weathergenr:::.is_int_scalar(2))
  testthat::expect_true(weathergenr:::.is_int_scalar(2.0))
  testthat::expect_false(weathergenr:::.is_int_scalar(2.5))
  testthat::expect_false(weathergenr:::.is_int_scalar(c(1, 2)))
  testthat::expect_false(weathergenr:::.is_int_scalar(NA_real_))
  testthat::expect_false(weathergenr:::.is_int_scalar(Inf))
  testthat::expect_false(weathergenr:::.is_int_scalar(TRUE))
})

testthat::test_that(".safe_cor handles insufficient and sufficient pairs", {
  x <- c(1, 2, NA, 4)
  y <- c(1, 2, 3, NA)
  out_low <- weathergenr:::.safe_cor(x, y, min_pairs = 3)
  testthat::expect_true(is.na(out_low["value"]))
  testthat::expect_equal(unname(out_low["n"]), 2)

  x2 <- c(1, 2, 3, 4)
  y2 <- c(1, 2, 3, 6)
  out_ok <- weathergenr:::.safe_cor(x2, y2, min_pairs = 3)
  testthat::expect_equal(unname(out_ok["n"]), 4)
  testthat::expect_equal(unname(out_ok["value"]), stats::cor(x2, y2))
})

testthat::test_that("format helpers handle NA and formatting", {
  testthat::expect_true(is.na(weathergenr:::.format_pct(NA_real_)))
  testthat::expect_true(is.na(weathergenr:::.format_num(NA_real_)))

  testthat::expect_identical(weathergenr:::.format_pct(1.234, digits = 1), "+1.2%")
  testthat::expect_identical(weathergenr:::.format_num(1.234, digits = 2), "1.23")
})

testthat::test_that(".log respects verbose and levels", {
  x <- 5
  testthat::expect_silent(weathergenr:::.log("hidden", verbose = FALSE))

  testthat::expect_message(
    weathergenr:::.log("value {x}", level = "info", verbose = TRUE),
    regexp = "value 5"
  )

  testthat::expect_warning(
    weathergenr:::.log("warn {x}", level = "warn", verbose = TRUE),
    regexp = "warn 5"
  )

  testthat::expect_error(
    weathergenr:::.log("error {x}", level = "error", verbose = TRUE),
    regexp = "error 5"
  )
})

testthat::test_that("format_elapsed returns unit-based strings", {
  out_secs <- weathergenr:::format_elapsed(Sys.time() - 30)
  out_mins <- weathergenr:::format_elapsed(Sys.time() - 90)
  out_hrs <- weathergenr:::format_elapsed(Sys.time() - 4000)

  testthat::expect_true(grepl("seconds", out_secs))
  testthat::expect_true(grepl("min", out_mins))
  testthat::expect_true(grepl("hr", out_hrs))
})
