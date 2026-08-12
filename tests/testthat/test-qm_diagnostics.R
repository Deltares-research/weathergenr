# ---- input validation --------------------------------------------------------

testthat::test_that("diagnose_precip_qm: mismatched series lengths are rejected", {
  testthat::expect_error(
    diagnose_precip_qm(precip_ref = c(1, 2, 3), precip_adj = c(1, 2)),
    "must have same length"
  )
})

testthat::test_that("diagnose_precip_qm: month and year must match the series length", {
  ref <- rep(1, 10)

  testthat::expect_error(
    diagnose_precip_qm(ref, ref, month = rep(1L, 9)),
    "'month' must have same length"
  )
  testthat::expect_error(
    diagnose_precip_qm(ref, ref, month = rep(1L, 10), year = rep(1L, 9)),
    "'year' must have same length"
  )
})

testthat::test_that("validate_quantile_mapping: wrapper rejects mismatched lengths before delegating", {
  testthat::expect_error(
    validate_quantile_mapping(precip_org = c(1, 2, 3), precip_adjusted = c(1, 2)),
    "'precip_org' and 'precip_adjusted' must have same length"
  )
})

# ---- structure ---------------------------------------------------------------

testthat::test_that("diagnose_precip_qm: optional blocks are NULL unless month/year are supplied", {
  set.seed(42)
  n_days <- 365L * 2L
  precip <- rgamma(n_days, shape = 1.2, scale = 5)
  precip[sample.int(n_days, size = 100)] <- 0

  bare <- diagnose_precip_qm(precip, precip)

  testthat::expect_s3_class(bare, "precip_qm_diagnostics")
  testthat::expect_null(bare$temporal)
  testthat::expect_null(bare$monthly)

  withmonth <- diagnose_precip_qm(
    precip, precip,
    month = rep(1:12, length.out = n_days)
  )

  # `monthly` needs only `month`; `temporal` needs both `month` and `year`.
  testthat::expect_false(is.null(withmonth$monthly))
  testthat::expect_null(withmonth$temporal)
})

testthat::test_that("diagnose_precip_qm: an identical series reports unit ratios and zero change", {
  set.seed(7)
  n_days <- 365L
  precip <- rgamma(n_days, shape = 1.5, scale = 4)
  precip[sample.int(n_days, size = 60)] <- 0

  out <- diagnose_precip_qm(precip, precip)

  moments <- out$moments
  testthat::expect_equal(moments$ratio, rep(1, nrow(moments)))
  testthat::expect_equal(moments$pct_change, rep(0, nrow(moments)))

  testthat::expect_equal(out$quantiles$ratio, rep(1, nrow(out$quantiles)))
  testthat::expect_equal(out$drydays$diff, rep(0, nrow(out$drydays)))
  testthat::expect_true(all(out$drydays$within_tolerance))
})

testthat::test_that("diagnose_precip_qm: intended ratios are recovered from the scenario factors", {
  set.seed(11)
  n_years <- 3L
  n_days <- 365L * n_years
  year <- rep(seq_len(n_years), each = 365L)
  month <- rep(1:12, length.out = n_days)

  precip <- rgamma(n_days, shape = 1.5, scale = 4)

  # A uniform +20% mean scenario applied exactly, so observed == intended.
  mean_factor <- matrix(1.2, nrow = n_years, ncol = 12)
  var_factor <- matrix(1.2^2, nrow = n_years, ncol = 12)

  out <- diagnose_precip_qm(
    precip_ref = precip,
    precip_adj = precip * 1.2,
    month = month,
    year = year,
    mean_factor = mean_factor,
    var_factor = var_factor
  )

  mean_row <- out$moments[out$moments$metric == "mean", ]
  testthat::expect_equal(mean_row$intended_ratio, 1.2)
  testthat::expect_equal(mean_row$ratio, 1.2)
  testthat::expect_equal(mean_row$ratio_error, 0)
  testthat::expect_true(mean_row$within_tolerance)
})

# ---- regression: NA propagation ----------------------------------------------
#
# Both cases below are guarded in the source with an explicit comment but had no
# test pinning the behaviour.

testthat::test_that("compute_spell_diagnostics: absent spells yield NA tolerance, not FALSE", {
  # Every day is wet, so there are no dry spells at all.
  all_wet <- rep(5, 100)

  out <- compute_spell_diagnostics(
    precip_ref = all_wet,
    precip_adj = all_wet,
    wet_thresh = 0.1
  )

  dry <- out[out$spell_type == "dry", ]
  testthat::expect_true(is.na(dry$mean_ratio))
  testthat::expect_true(is.na(dry$within_tolerance))
  testthat::expect_false(isFALSE(dry$within_tolerance))

  # The wet row is well defined and unchanged, so it must still pass.
  wet <- out[out$spell_type == "wet", ]
  testthat::expect_equal(wet$mean_ratio, 1)
  testthat::expect_true(wet$within_tolerance)
})

testthat::test_that("compute_dryday_diagnostics: a zero original count does not error on the ratio", {
  # No dry days in either series, so `dry_days` has original == 0.
  all_wet <- rep(5, 50)

  out <- testthat::expect_no_warning(
    compute_dryday_diagnostics(all_wet, all_wet, wet_thresh = 0.1)
  )

  dry <- out[out$category == "dry_days", ]
  testthat::expect_equal(dry$original, 0)
  testthat::expect_equal(dry$diff, 0)
  testthat::expect_true(dry$within_tolerance)
  testthat::expect_true(is.nan(dry$ratio) || is.na(dry$ratio))
})

# ---- print / summary snapshots -----------------------------------------------
#
# Deterministic hand-built inputs, not RNG: these snapshots pin exact formatted
# numbers, which must not shift with an R version's RNG stream.

testthat::test_that("print and summary methods are stable", {
  # A distinct level per month, so the monthly correlation is well defined
  # (a flat month-to-month profile gives a zero-variance `cor()` input).
  base <- c(0, 1, 2, 4, 8, 0, 3, 6, 0, 5)
  seasonal <- c(1, 1.5, 2, 2.5, 3, 3.5, 4, 3.5, 3, 2.5, 2, 1.5)

  ref <- as.numeric(outer(base, seasonal))

  # Deliberately 1.2, not 1.25. The snapshot prints the achieved variance ratio
  # at three decimals, and 1.25^2 is exactly 1.5625 -- an exact rounding
  # midpoint, which C libraries break differently: macOS rounds half to even
  # ("1.562"), Windows and glibc round half away from zero ("1.563"). That made
  # this snapshot fail on macOS only. 1.2^2 = 1.44 is not a midpoint, so the
  # printed digits are stable everywhere.
  adj <- ref * 1.2

  n <- length(ref)
  month <- rep(1:12, each = length(base))
  year <- rep(1L, n)

  out <- diagnose_precip_qm(
    precip_ref = ref,
    precip_adj = adj,
    month = month,
    year = year,
    mean_factor = matrix(1.2, nrow = 1, ncol = 12),
    var_factor = matrix(1.2^2, nrow = 1, ncol = 12)
  )

  testthat::expect_snapshot(print(out))
  testthat::expect_snapshot(summary(out))
})

testthat::test_that("print and summary return their input invisibly", {
  ref <- rep(c(0, 1, 2, 4, 8), times = 20)
  out <- diagnose_precip_qm(ref, ref * 1.1)

  utils::capture.output({
    testthat::expect_invisible(print(out))
    testthat::expect_invisible(summary(out))
    testthat::expect_identical(print(out), out)
  })
})
