

# =============================================================================
# Unit Tests for wavelet_spectral_analysis
# =============================================================================

library(testthat)
library(weathergenr)

# =============================================================================
# BASIC STRUCTURE AND INPUT VALIDATION
# =============================================================================

test_that("wavelet_spectral_analysis returns expected structure", {
  set.seed(123)
  x <- rnorm(64)

  res <- wavelet_spectral_analysis(x)

  expect_type(res, "list")

  required <- c(
    "gws", "gws_signif", "gws_period",
    "signif_periods", "comps", "coi",
    "has_significance", "significance_status"
  )

  expect_true(all(required %in% names(res)))
  expect_true(is.numeric(res$gws))
  expect_true(is.numeric(res$gws_signif))
  expect_true(is.numeric(res$gws_period))
  expect_true(is.integer(res$signif_periods))
  expect_true(is.matrix(res$comps))
})

test_that("wavelet_spectral_analysis rejects NA input", {
  x <- rnorm(64)
  x[10] <- NA

  expect_error(
    wavelet_spectral_analysis(x),
    "missing values"
  )
})

test_that("wavelet_spectral_analysis is deterministic", {
  set.seed(999)
  x <- rnorm(64)

  res1 <- wavelet_spectral_analysis(x)
  res2 <- wavelet_spectral_analysis(x)

  expect_equal(res1$gws, res2$gws)
  expect_equal(res1$gws_signif, res2$gws_signif)
  expect_equal(res1$signif_periods, res2$signif_periods)
})

# =============================================================================
# RECONSTRUCTION / CLOSURE TESTS (CORE CONTRACT)
# =============================================================================

test_that("comps exactly reconstruct the original signal (multiple cases)", {

  cases <- list(
    pure_sine = {
      n <- 128
      2 * sin(2 * pi * (1:n) / 8)
    },
    mixed_signal = {
      n <- 256
      3 * sin(2 * pi * (1:n) / 8) +
        1.5 * sin(2 * pi * (1:n) / 16) +
        rnorm(n, sd = 0.5)
    },
    ar1_process = {
      as.numeric(arima.sim(n = 128, model = list(ar = 0.7)))
    }
  )

  for (signal in cases) {
    res <- wavelet_spectral_analysis(signal, warn_neff = FALSE)

    recon <- rowSums(res$comps)
    err <- max(abs(signal - recon))

    expect_lt(
      err,
      1e-6,
      label = "Reconstruction closure failed"
    )
  }
})

test_that("comps exists and equals NOISE-only when no significance is found", {
  set.seed(1)
  x <- rnorm(64)

  res <- wavelet_spectral_analysis(
    x,
    noise.type = "white",
    signif.level = 0.99
  )

  if (!res$has_significance) {
    expect_equal(colnames(res$comps), "NOISE")
    expect_equal(res$comps[, 1], x)
  }
})

# =============================================================================
# SIGNAL DETECTION (SOFT ASSERTIONS)
# =============================================================================

test_that("detects known periodic signal", {
  set.seed(123)
  n <- 128
  t <- 1:n
  true_period <- 8

  x <- sin(2 * pi * t / true_period) + rnorm(n, sd = 0.3)

  res <- wavelet_spectral_analysis(x)

  expect_true(length(res$signif_periods) > 0)

  detected <- res$gws_period[res$signif_periods]

  expect_true(
    any(abs(detected - true_period) < 0.25 * true_period),
    label = "Known periodicity not detected"
  )
})

test_that("white noise does not systematically produce long-period significance", {
  set.seed(42)
  x <- rnorm(128)

  res <- wavelet_spectral_analysis(x, noise.type = "white")

  if (length(res$signif_periods) > 0) {
    periods <- res$gws_period[res$signif_periods]
    expect_true(
      any(periods < length(x) / 2),
      label = "Unexpected dominance of very long periods in white noise"
    )
  }
})

# =============================================================================
# COI / POINTWISE SIGNIFICANCE SANITY
# =============================================================================

test_that("COI masking removes pointwise significance outside COI", {
  set.seed(42)
  n <- 64
  t <- seq_len(n)

  x <- sin(2 * pi * t / n)

  res <- wavelet_spectral_analysis(x, noise.type = "white")

  expect_true("sigm_coi" %in% names(res))

  # Outside COI must be NA
  expect_true(all(is.na(res$sigm_coi[!outer(res$gws_period, res$coi, "<=")])))
})

# =============================================================================
# LAG-1 BOOTSTRAP DIAGNOSTICS
# =============================================================================

test_that("lag1 bootstrap CI is returned when requested", {
  set.seed(101)
  x <- as.numeric(arima.sim(n = 128, model = list(ar = 0.7)))

  res <- wavelet_spectral_analysis(
    x,
    lag1_ci = TRUE,
    lag1_ci_level = 0.90,
    lag1_boot_n = 200,
    seed = 999,
    warn_neff = FALSE
  )

  expect_true("diagnostics" %in% names(res))
  ci <- res$diagnostics$lag1_ci

  expect_true(is.list(ci))
  expect_true(all(c("lag1_hat", "lower", "upper") %in% names(ci)))
  expect_true(ci$lower <= ci$lag1_hat && ci$lag1_hat <= ci$upper)
})

# =============================================================================
# Neff WARNING BEHAVIOR
# =============================================================================

test_that("warning is emitted when Neff is small for most scales", {
  set.seed(202)
  x <- as.numeric(arima.sim(n = 16, model = list(ar = 0.6)))

  expect_warning(
    wavelet_spectral_analysis(
      x,
      warn_neff = TRUE,
      neff_warn_min = 10,
      neff_warn_frac = 0.5
    ),
    regexp = "neff"
  )
})

# =============================================================================
# DIAGNOSTICS OUTPUT
# =============================================================================

test_that("diagnostics contain required fields", {
  set.seed(123)
  x <- rnorm(64)

  res <- wavelet_spectral_analysis(x)

  diag <- res$diagnostics

  required <- c(
    "lag1", "variance", "n_original",
    "scale", "fourier_factor",
    "n_coi", "gws_neff"
  )

  expect_true(all(required %in% names(diag)))
  expect_true(diag$variance > 0)
  expect_equal(diag$n_original, 64)
})
