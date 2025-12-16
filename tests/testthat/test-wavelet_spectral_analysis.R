
library(testthat)

test_that("wavelet_spectral_analysis returns expected list structure", {

  set.seed(123)
  x <- rnorm(64)

  res <- wavelet_spectral_analysis(x, plot = FALSE)

  expect_type(res, "list")
  expect_true(all(c("GWS", "GWS_signif", "GWS_period", "signif_periods") %in% names(res)))

  expect_true(is.numeric(res$GWS))
  expect_true(is.numeric(res$GWS_signif))
  expect_true(is.numeric(res$GWS_period))
  expect_true(is.integer(res$signif_periods) || is.numeric(res$signif_periods))
})

test_that("wavelet_spectral_analysis rejects non-numeric input", {
  expect_error(
    wavelet_spectral_analysis(c("a", "b", "c")),
    "'variable' must be numeric"
  )
})

test_that("wavelet_spectral_analysis rejects NA values", {
  x <- rnorm(32)
  x[10] <- NA
  expect_error(
    wavelet_spectral_analysis(x),
    "missing values"
  )
})

test_that("wavelet_spectral_analysis rejects short time series", {
  expect_error(
    wavelet_spectral_analysis(rnorm(8)),
    "'variable' must be longer than 8 values"
  )
})

test_that("wavelet_spectral_analysis validates arguments", {
  x <- rnorm(64)

  expect_error(
    wavelet_spectral_analysis(x, signif.level = 1.5),
    "'signif.level' must be a finite numeric value in \\(0, 1\\)"
  )
})


test_that("white noise does not produce strong low-frequency signal", {

  set.seed(42)
  x <- rnorm(128)

  res <- wavelet_spectral_analysis(x, noise.type = "white")

  if (length(res$signif_periods) > 0) {
    periods <- res$GWS_period[res$signif_periods]
    expect_true(all(periods < length(x) / 2))
  }
})

test_that("wavelet_spectral_analysis detects known periodic signal", {

  set.seed(123)
  n <- 128
  t <- 1:n
  true_period <- 8

  x <- sin(2 * pi * t / true_period) + rnorm(n, sd = 0.3)

  res <- wavelet_spectral_analysis(x, signif.level = 0.90)

  expect_true(length(res$signif_periods) > 0)

  detected_periods <- res$GWS_period[res$signif_periods]

  # Allow tolerance of +/- 20%
  expect_true(any(abs(detected_periods - true_period) < 0.2 * true_period))
})


test_that("wavelet components reconstruct original signal", {

  set.seed(123)
  n <- 64
  t <- 1:n
  x <- sin(2 * pi * t / 8) + rnorm(n, sd = 0.2)

  res <- wavelet_spectral_analysis(x)

  if (!is.null(res$COMPS)) {

    comps <- as.matrix(res$COMPS[, grepl("^Component_", names(res$COMPS))])
    noise <- res$COMPS$NOISE

    recon <- rowSums(comps) + noise

    # Reconstruction should be close (not exact)
    expect_lt(
      mean((x - recon)^2),
      var(x)
    )
  }
})

test_that("wavelet_spectral_analysis is deterministic", {

  set.seed(999)
  x <- rnorm(64)

  res1 <- wavelet_spectral_analysis(x)
  res2 <- wavelet_spectral_analysis(x)

  expect_equal(res1$GWS, res2$GWS)
  expect_equal(res1$GWS_signif, res2$GWS_signif)
  expect_equal(res1$signif_periods, res2$signif_periods)
})

test_that("wavelet_spectral_analysis plotting does not error", {

  skip_on_cran()

  set.seed(123)
  x <- rnorm(64)

  tmp <- tempdir()

  expect_silent(
    wavelet_spectral_analysis(
      x,
      plot = TRUE,
      output.path = tmp
    )
  )
})
test_that("wavelet_spectral_analysis enforces COI masking for low-frequency power", {

  set.seed(42)

  n <- 64
  t <- seq_len(n)

  # Strong low-frequency + edge-dominated signal
  # This WILL trigger false low-frequency detection if COI is not enforced
  x <- sin(2 * pi * t / n) + 0.01 * t

  res <- wavelet_spectral_analysis(
    variable = x,
    signif.level = 0.90,
    noise.type = "white",
    plot = FALSE
  )

  # Basic sanity
  expect_type(res, "list")
  expect_true("signif_periods" %in% names(res))

  # Critical assertion:
  # No detected period should be close to the full series length
  # (COI masking should remove this)
  if (length(res$signif_periods) > 0) {
    expect_true(
      all(res$signif_periods < n / 2),
      info = "COI masking failed: detected periods near series length"
    )
  } else {
    succeed("No significant periods detected (acceptable under COI masking)")
  }
})

