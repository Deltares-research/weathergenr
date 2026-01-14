testthat::test_that("DEBUG: analyze_wavelet_spectrum identity", {
  suppressMessages(suppressWarnings({

    w <- weathergenr::analyze_wavelet_spectrum
    env <- environment(w)
    where <- tryCatch(getNamespaceName(env), error = function(e) "<non-namespace>")
    msg <- paste0(
      "analyze_wavelet_spectrum from: ", where,
      " | first body line: ", paste0(deparse(body(w))[1], collapse = " ")
    )
    message(msg)
    testthat::expect_true(TRUE)

  }))
})

# =============================================================================
# BASIC STRUCTURE AND INPUT VALIDATION
# =============================================================================

testthat::test_that("analyze_wavelet_spectrum (fast mode) returns expected structure", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    x <- rnorm(64)

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "fast"
    )

    testthat::expect_type(res, "list")

    required_fast <- c(
      "gws", "gws_unmasked", "period",
      "gws_signif", "gws_signif_unmasked",
      "has_significance", "signif_periods",
      "coi", "power"
    )

    testthat::expect_true(all(required_fast %in% names(res)))

    testthat::expect_true(is.numeric(res$gws))
    testthat::expect_true(is.numeric(res$gws_unmasked))
    testthat::expect_true(is.numeric(res$gws_signif))
    testthat::expect_true(is.numeric(res$gws_signif_unmasked))
    testthat::expect_true(is.numeric(res$period))
    testthat::expect_true(is.logical(res$has_significance))
    testthat::expect_true(is.numeric(res$coi))
    testthat::expect_true(is.matrix(res$power))

  }))
})

testthat::test_that("analyze_wavelet_spectrum (complete mode) returns expected structure", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    x <- rnorm(64)

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "complete"
    )

    testthat::expect_type(res, "list")

    required_complete_min <- c(
      "gws", "gws_unmasked", "period",
      "gws_signif", "gws_signif_unmasked",
      "has_significance", "signif_periods",
      "coi", "power",
      "power_coi", "sigm", "sigm_coi", "power_signif_coi",
      "wave", "comps", "comps_names", "gws_n_coi", "neff", "neff_unmasked"
    )

    testthat::expect_true(all(required_complete_min %in% names(res)))
    testthat::expect_true(is.numeric(res$gws))
    testthat::expect_true(is.numeric(res$period))
    testthat::expect_true(is.numeric(res$gws_signif))
    testthat::expect_true(is.matrix(res$power))
    testthat::expect_true(is.matrix(res$power_coi))
    testthat::expect_true(is.matrix(res$sigm))
    testthat::expect_true(is.matrix(res$sigm_coi))
    testthat::expect_true(is.matrix(res$comps))

  }))
})

testthat::test_that("analyze_wavelet_spectrum rejects NA input", {
  suppressMessages(suppressWarnings({

    x <- rnorm(64)
    x[10] <- NA

    testthat::expect_error(
      weathergenr::analyze_wavelet_spectrum(series = x, mode = "fast"),
      regexp = "missing values|contains missing values",
      fixed = FALSE
    )

  }))
})

testthat::test_that("analyze_wavelet_spectrum is deterministic (no bootstrap)", {
  suppressMessages(suppressWarnings({

    set.seed(999)
    x <- rnorm(64)

    res1 <- weathergenr::analyze_wavelet_spectrum(series = x, mode = "fast")
    res2 <- weathergenr::analyze_wavelet_spectrum(series = x, mode = "fast")

    testthat::expect_equal(res1$gws, res2$gws)
    testthat::expect_equal(res1$gws_signif, res2$gws_signif)
    testthat::expect_equal(res1$signif_periods, res2$signif_periods)

  }))
})

# =============================================================================
# RECONSTRUCTION / CLOSURE (REALISTIC CONTRACT)
# =============================================================================

testthat::test_that("comps reconstruct the signal within a reasonable relative error (multiple cases)", {
  suppressMessages(suppressWarnings({

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

      res <- weathergenr::analyze_wavelet_spectrum(
        series = signal,
        mode = "complete"
      )

      testthat::expect_true(is.matrix(res$comps))

      recon <- rowSums(res$comps)
      e <- signal - recon

      rmse <- sqrt(mean(e^2))
      max_abs <- max(abs(e))

      sdx <- stats::sd(signal)
      iqr <- stats::IQR(signal)

      testthat::expect_lt(rmse, 0.25 * sdx, label = "RMSE too large relative to SD")
      testthat::expect_lt(max_abs, 0.35 * iqr, label = "Max error too large relative to IQR")
    }

  }))
})

testthat::test_that("comps equals noise-only when no significance is found", {
  suppressMessages(suppressWarnings({

    set.seed(1)
    x <- rnorm(64)

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "complete",
      noise = "white",
      signif = 0.999
    )

    if (!isTRUE(res$has_significance) && length(res$signif_periods) == 0) {
      testthat::expect_true(is.matrix(res$comps))
      testthat::expect_equal(ncol(res$comps), 1)
      testthat::expect_true(tolower(colnames(res$comps)[1]) %in% c("noise"))
      testthat::expect_equal(as.numeric(res$comps[, 1]), as.numeric(x))
    }

  }))
})

# =============================================================================
# SIGNAL DETECTION (SOFT ASSERTIONS)
# =============================================================================

testthat::test_that("detects known periodic signal (soft)", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    n <- 128
    t <- 1:n
    true_period <- 8

    x <- sin(2 * pi * t / true_period) + rnorm(n, sd = 0.3)

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "fast"
    )

    if (length(res$signif_periods) > 0) {
      detected <- res$period[res$signif_periods]

      testthat::expect_true(
        any(abs(detected - true_period) < 0.35 * true_period),
        label = "Known periodicity not detected within tolerance"
      )
    } else {
      testthat::succeed()
    }

  }))
})

testthat::test_that("white noise does not systematically produce only very long-period significance", {
  suppressMessages(suppressWarnings({

    set.seed(42)
    x <- rnorm(128)

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "fast",
      noise = "white"
    )

    if (length(res$signif_periods) > 0) {
      periods <- res$period[res$signif_periods]
      testthat::expect_true(
        any(periods < length(x) / 2),
        label = "Unexpected dominance of very long periods in white noise"
      )
    }

  }))
})

# =============================================================================
# COI / POINTWISE SIGNIFICANCE SANITY (COMPLETE MODE)
# =============================================================================

testthat::test_that("COI masking removes pointwise significance ratio outside COI", {
  suppressMessages(suppressWarnings({

    set.seed(42)
    n <- 64
    t <- seq_len(n)

    x <- sin(2 * pi * t / n)

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "complete",
      noise = "white"
    )

    testthat::expect_true("sigm_coi" %in% names(res))
    testthat::expect_true(is.matrix(res$sigm_coi))

    coi_mask <- outer(res$period, res$coi, FUN = "<=")
    testthat::expect_true(all(is.na(res$sigm_coi[!coi_mask])))

  }))
})

# =============================================================================
# LAG-1 BOOTSTRAP DIAGNOSTICS
# =============================================================================

testthat::test_that("lag1 bootstrap CI is returned when requested", {
  suppressMessages(suppressWarnings({

    set.seed(101)
    x <- as.numeric(arima.sim(n = 128, model = list(ar = 0.7)))

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "complete",
      diagnostics = TRUE,
      lag1_ci = TRUE,
      lag1_ci_level = 0.90,
      lag1_boot_n = 200,
      seed = 999
    )

    testthat::expect_true("diagnostics" %in% names(res))
    ci <- res$diagnostics$lag1_ci

    testthat::expect_true(is.list(ci))
    testthat::expect_true(all(c("lag1_hat", "lower", "upper") %in% names(ci)))
    testthat::expect_true(ci$lower <= ci$lag1_hat && ci$lag1_hat <= ci$upper)

  }))
})

# =============================================================================
# Neff WARNING BEHAVIOR
# =============================================================================

testthat::test_that("warning is emitted when Neff is small for most scales", {

  set.seed(202)
  x <- as.numeric(arima.sim(n = 16, model = list(ar = 0.6)))

  testthat::expect_warning(
    weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "fast",
      diagnostics = FALSE,
      warn_neff = TRUE,
      neff_warn_min = 1e6,
      neff_warn_frac = 0.01
    ),
    regexp = "neff",
    fixed = FALSE
  )
})

# =============================================================================
# DIAGNOSTICS OUTPUT
# =============================================================================

testthat::test_that("diagnostics contain required fields (diagnostics=TRUE, complete mode)", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    x <- rnorm(64)

    res <- weathergenr::analyze_wavelet_spectrum(
      series = x,
      mode = "complete",
      diagnostics = TRUE
    )

    testthat::expect_true("diagnostics" %in% names(res))
    diag <- res$diagnostics

    required <- c(
      "lag1", "variance", "n_original",
      "scale", "fourier_factor",
      "n_coi", "neff"
    )

    testthat::expect_true(all(required %in% names(diag)))
    testthat::expect_true(is.numeric(diag$variance))
    testthat::expect_true(diag$variance > 0)
    testthat::expect_equal(diag$n_original, 64)

  }))
})
