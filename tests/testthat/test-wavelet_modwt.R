# Functions tested (relative paths):
# - R/wavelet_modwt.R: modwt_decompose(), modwt_reconstruct(), modwt_mra(),
#   analyze_wavelet_additive(), print.wavelet_additive()

testthat::test_that("modwt_decompose validates inputs and returns structure", {
  series <- sin(seq(0, 2 * pi, length.out = 32))

  out <- modwt_decompose(series, filter = "la8", n_levels = 2)
  testthat::expect_s3_class(out, "modwt_result")
  testthat::expect_true(all(c("coefficients", "filter", "n_levels", "n",
                              "boundary", "filter_length") %in% names(out)))
  testthat::expect_equal(out$n, length(series))
  testthat::expect_true(is.list(out$coefficients))
  testthat::expect_equal(length(out$coefficients), out$n_levels + 1L)
  testthat::expect_equal(
    names(out$coefficients),
    c(paste0("d", seq_len(out$n_levels)), paste0("s", out$n_levels))
  )

  testthat::expect_warning(
    modwt_decompose(series, n_levels = 99),
    "exceeds recommended"
  )

  testthat::expect_error(modwt_decompose(1:7), "at least 8")
  testthat::expect_error(modwt_decompose(c(1, NA, 2)), "missing values")
  testthat::expect_error(modwt_decompose("a"), "numeric vector")
  testthat::expect_error(modwt_decompose(1:8, boundary = "reflect"), "periodic")
})

testthat::test_that("modwt_reconstruct handles modwt_result and raw lists", {
  series <- sin(seq(0, 2 * pi, length.out = 32))
  wt <- modwt_decompose(series, n_levels = 2)

  recon <- modwt_reconstruct(wt)
  testthat::expect_equal(length(recon), length(series))
  testthat::expect_lt(max(abs(recon - series)), 1e-6)

  testthat::expect_error(modwt_reconstruct(list(a = 1)), "requires W")
  testthat::expect_error(modwt_reconstruct(1), "modwt_result")
})

testthat::test_that("modwt_mra returns additive components and diagnostics", {
  series <- sin(seq(0, 2 * pi, length.out = 32))

  out <- modwt_mra(series, n_levels = 2, include_smooth = TRUE)
  testthat::expect_true(is.matrix(out$components))
  testthat::expect_equal(dim(out$components), c(length(series), out$n_levels + 1L))
  testthat::expect_equal(colnames(out$components), c("D1", "D2", "S2"))
  testthat::expect_equal(length(out$periods), ncol(out$components))
  testthat::expect_equal(length(out$variance), ncol(out$components))
  testthat::expect_equal(length(out$variance_fraction), ncol(out$components))
  testthat::expect_true(isTRUE(out$is_additive))
  testthat::expect_lt(out$reconstruction_error, 1e-6)

  no_smooth <- modwt_mra(series, n_levels = 2, include_smooth = FALSE)
  testthat::expect_equal(ncol(no_smooth$components), no_smooth$n_levels)
  testthat::expect_equal(colnames(no_smooth$components), c("D1", "D2"))

  testthat::expect_error(modwt_mra(series, max_period_frac = 0), "in \\(0, 1\\]")
})

testthat::test_that("modwt_mra handles constant series deterministically", {
  series <- rep(5, 32)
  out <- modwt_mra(series, n_levels = 2)

  testthat::expect_equal(out$total_variance, 0)
  testthat::expect_true(all(out$variance == 0))
  testthat::expect_lt(out$reconstruction_error, 1e-8)
})

testthat::test_that("analyze_wavelet_additive returns expected structure", {
  series <- sin(seq(0, 2 * pi, length.out = 32))

  out <- analyze_wavelet_additive(
    series,
    noise = "white",
    cwt_mode = "fast",
    diagnostics = TRUE
  )

  testthat::expect_s3_class(out, "wavelet_additive")
  testthat::expect_true(is.matrix(out$components))
  testthat::expect_equal(nrow(out$components), length(series))
  testthat::expect_equal(length(out$periods), ncol(out$components))
  testthat::expect_equal(length(out$variance), ncol(out$components))
  testthat::expect_true(is.list(out$diagnostics))

  testthat::expect_error(
    analyze_wavelet_additive(series[1:10]),
    "at least 16"
  )
  testthat::expect_error(
    analyze_wavelet_additive(c(1, NA, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)),
    "missing values"
  )
})

testthat::test_that("print.wavelet_additive prints a compact summary", {
  series <- sin(seq(0, 2 * pi, length.out = 32))
  out <- analyze_wavelet_additive(series, noise = "white", cwt_mode = "fast")

  testthat::expect_output(print(out), "Additive Wavelet Analysis")
})
