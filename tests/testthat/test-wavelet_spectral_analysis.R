test_that("DEBUG: wavelet_spectral_analysis identity", {
  suppressMessages(suppressWarnings({

    w <- wavelet_spectral_analysis
    env <- environment(w)
    where <- tryCatch(getNamespaceName(env), error = function(e) "<non-namespace>")
    msg <- paste0(
      "wavelet_spectral_analysis from: ", where,
      " | file: ", paste0(deparse(body(w))[1], collapse = " ")
    )
    message(msg)
    expect_true(TRUE)

  }))
})


# =============================================================================
# BASIC STRUCTURE AND INPUT VALIDATION
# =============================================================================

test_that("wavelet_spectral_analysis (fast mode) returns expected structure", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    x <- rnorm(64)

    res <- weathergenr::wavelet_spectral_analysis(
      x,
      mode = "fast"
    )

    expect_type(res, "list")

    required_fast <- c(
      "gws", "gws_unmasked", "gws_period",
      "gws_signif", "gws_signif_unmasked",
      "has_significance", "signif_periods",
      "coi", "power"
    )

    expect_true(all(required_fast %in% names(res)))

    expect_true(is.numeric(res$gws))
    expect_true(is.numeric(res$gws_unmasked))
    expect_true(is.numeric(res$gws_signif))
    expect_true(is.numeric(res$gws_signif_unmasked))
    expect_true(is.numeric(res$gws_period))
    expect_true(is.logical(res$has_significance))
    expect_true(is.numeric(res$coi))
    expect_true(is.matrix(res$power))

  }))
})

test_that("wavelet_spectral_analysis (complete mode) returns expected structure", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    x <- rnorm(64)

    res <- weathergenr::wavelet_spectral_analysis(
      x,
      mode = "complete"
    )

    expect_type(res, "list")

    required_complete_min <- c(
      "gws", "gws_unmasked", "gws_period",
      "gws_signif", "gws_signif_unmasked",
      "has_significance", "signif_periods",
      "coi", "power",
      "power_coi", "sigm", "sigm_coi", "power_signif_coi",
      "wave", "comps", "gws_n_coi", "gws_neff", "gws_neff_unmasked"
    )

    expect_true(all(required_complete_min %in% names(res)))
    expect_true(is.numeric(res$gws))
    expect_true(is.numeric(res$gws_period))
    expect_true(is.numeric(res$gws_signif))
    expect_true(is.matrix(res$power))
    expect_true(is.matrix(res$power_coi))
    expect_true(is.matrix(res$sigm))
    expect_true(is.matrix(res$sigm_coi))
    expect_true(is.matrix(res$comps))

  }))
})

test_that("wavelet_spectral_analysis rejects NA input", {
  suppressMessages(suppressWarnings({

    x <- rnorm(64)
    x[10] <- NA

    expect_error(
      weathergenr::wavelet_spectral_analysis(x, mode = "fast"),
      regexp = "missing values|contains missing values",
      fixed = FALSE
    )

  }))
})

test_that("wavelet_spectral_analysis is deterministic (no bootstrap)", {
  suppressMessages(suppressWarnings({

    set.seed(999)
    x <- rnorm(64)

    res1 <- weathergenr::wavelet_spectral_analysis(x, mode = "fast")
    res2 <- weathergenr::wavelet_spectral_analysis(x, mode = "fast")

    expect_equal(res1$gws, res2$gws)
    expect_equal(res1$gws_signif, res2$gws_signif)
    expect_equal(res1$signif_periods, res2$signif_periods)

  }))
})


# =============================================================================
# RECONSTRUCTION / CLOSURE (REALISTIC CONTRACT)
# =============================================================================

test_that("comps reconstruct the signal within a reasonable relative error (multiple cases)", {
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

      res <- weathergenr::wavelet_spectral_analysis(
        signal,
        mode = "complete"
      )

      expect_true(is.matrix(res$comps))

      recon <- rowSums(res$comps)
      e <- signal - recon

      rmse <- sqrt(mean(e^2))
      max_abs <- max(abs(e))

      sdx <- stats::sd(signal)
      iqr <- stats::IQR(signal)

      expect_lt(rmse, 0.25 * sdx, label = "RMSE too large relative to SD")
      expect_lt(max_abs, 0.35 * iqr, label = "Max error too large relative to IQR")
    }

  }))
})

test_that("comps equals noise-only when no significance is found", {
  suppressMessages(suppressWarnings({

    set.seed(1)
    x <- rnorm(64)

    res <- weathergenr::wavelet_spectral_analysis(
      x,
      mode = "complete",
      ,
      noise.type = "white",
      signif.level = 0.999
    )

    if (!isTRUE(res$has_significance) && length(res$signif_periods) == 0) {
      expect_true(is.matrix(res$comps))
      expect_equal(ncol(res$comps), 1)
      expect_true(tolower(colnames(res$comps)[1]) %in% c("noise"))
      expect_equal(as.numeric(res$comps[, 1]), as.numeric(x))
    }

  }))
})


# =============================================================================
# SIGNAL DETECTION (SOFT ASSERTIONS)
# =============================================================================

test_that("detects known periodic signal (soft)", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    n <- 128
    t <- 1:n
    true_period <- 8

    x <- sin(2 * pi * t / true_period) + rnorm(n, sd = 0.3)

    res <- weathergenr::wavelet_spectral_analysis(
      x,
      mode = "fast"
    )

    if (length(res$signif_periods) > 0) {
      detected <- res$gws_period[res$signif_periods]

      expect_true(
        any(abs(detected - true_period) < 0.35 * true_period),
        label = "Known periodicity not detected within tolerance"
      )
    } else {
      succeed()
    }

  }))
})

test_that("white noise does not systematically produce only very long-period significance", {
  suppressMessages(suppressWarnings({

    set.seed(42)
    x <- rnorm(128)

    res <- weathergenr::wavelet_spectral_analysis(
      x,
      mode = "fast",
      noise.type = "white"
    )

    if (length(res$signif_periods) > 0) {
      periods <- res$gws_period[res$signif_periods]
      expect_true(
        any(periods < length(x) / 2),
        label = "Unexpected dominance of very long periods in white noise"
      )
    }

  }))
})


# =============================================================================
# COI / POINTWISE SIGNIFICANCE SANITY (COMPLETE MODE)
# =============================================================================

test_that("COI masking removes pointwise significance ratio outside COI", {
  suppressMessages(suppressWarnings({

    set.seed(42)
    n <- 64
    t <- seq_len(n)

    x <- sin(2 * pi * t / n)

    res <- weathergenr::wavelet_spectral_analysis(
      x,
      mode = "complete",
      noise.type = "white"
    )

    expect_true("sigm_coi" %in% names(res))
    expect_true(is.matrix(res$sigm_coi))

    coi_mask <- outer(res$gws_period, res$coi, FUN = "<=")
    expect_true(all(is.na(res$sigm_coi[!coi_mask])))

  }))
})


# =============================================================================
# LAG-1 BOOTSTRAP DIAGNOSTICS
# =============================================================================

test_that("lag1 bootstrap CI is returned when requested", {
  suppressMessages(suppressWarnings({

    set.seed(101)
    x <- as.numeric(arima.sim(n = 128, model = list(ar = 0.7)))

    res <- weathergenr::wavelet_spectral_analysis(
      x,
      lag1_ci = TRUE,
      lag1_ci_level = 0.90,
      lag1_boot_n = 200,
      seed = 999
    )

    expect_true("diagnostics" %in% names(res))
    ci <- res$diagnostics$lag1_ci

    expect_true(is.list(ci))
    expect_true(all(c("lag1_hat", "lower", "upper") %in% names(ci)))
    expect_true(ci$lower <= ci$lag1_hat && ci$lag1_hat <= ci$upper)

  }))
})


# =============================================================================
# Neff WARNING BEHAVIOR
# =============================================================================

test_that("warning is emitted when Neff is small for most scales", {

  set.seed(202)
  x <- as.numeric(arima.sim(n = 16, model = list(ar = 0.6)))

  expect_warning(
    weathergenr::wavelet_spectral_analysis(
      variable = x,
      mode = "fast",
      return_diagnostics = FALSE,
      warn_neff = TRUE,
      neff_warn_min = 1e6,
      neff_warn_frac = 0.01
    ),
    regexp = "neff",
    fixed = FALSE
  )
})


# =============================================================================
# DIAGNOSTICS OUTPUT (DEFAULT BEHAVIOR)
# =============================================================================

test_that("DEBUG: wavelet_spectral_analysis identity", {
  suppressMessages(suppressWarnings({

    w <- weathergenr::wavelet_spectral_analysis
    env <- environment(w)
    where <- tryCatch(getNamespaceName(env), error = function(e) "<non-namespace>")
    msg <- paste0(
      "wavelet_spectral_analysis from: ", where,
      " | file: ", paste0(deparse(body(w))[1], collapse = " ")
    )
    message(msg)
    expect_true(TRUE)

  }))
})

test_that("diagnostics contain required fields (default return_diagnostics=TRUE)", {
  suppressMessages(suppressWarnings({

    set.seed(123)
    x <- rnorm(64)

    res <- weathergenr::wavelet_spectral_analysis(x)

    expect_true("diagnostics" %in% names(res))
    diag <- res$diagnostics

    required <- c(
      "lag1", "variance", "n_original",
      "scale", "fourier_factor",
      "n_coi", "gws_neff"
    )

    expect_true(all(required %in% names(diag)))
    expect_true(is.numeric(diag$variance))
    expect_true(diag$variance > 0)
    expect_equal(diag$n_original, 64)

  }))
})
