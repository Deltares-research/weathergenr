# test-filter_warm_realizations.R

library(testthat)
library(weathergenr)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

make_inputs <- function(n_years = 30, n_real = 10, n_period = 12) {
  set.seed(1)

  series.obs <- rnorm(n_years, mean = 100, sd = 20)
  series.sim <- replicate(n_real, rnorm(n_years, mean = 100, sd = 20))
  series.sim <- matrix(series.sim, nrow = n_years, ncol = n_real)

  power.period <- seq_len(n_period)
  power.obs <- runif(n_period, 0.5, 2.0)
  power.signif <- rep(1.0, n_period)

  # power.sim must be [n_period x n_real]
  power.sim <- matrix(runif(n_period * n_real, 0.2, 2.5), nrow = n_period, ncol = n_real)

  list(
    series.obs = series.obs,
    series.sim = series.sim,
    power.obs = power.obs,
    power.sim = power.sim,
    power.period = power.period,
    power.signif = power.signif
  )
}

# Capture all warnings produced by expr; return list(value=..., warnings=character())
capture_warnings <- function(expr) {
  warns <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = warns)
}

# -----------------------------------------------------------------------------
# 1) Basic structure + determinism
# -----------------------------------------------------------------------------

test_that("filter_warm_simulations returns expected structure", {
  inp <- make_inputs(n_years = 30, n_real = 12, n_period = 10)

  out <- filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    power.obs = inp$power.obs,
    power.sim = inp$power.sim,
    power.period = inp$power.period,
    power.signif = inp$power.signif,
    sample.num = 5,
    seed = 123,
    save.plots = FALSE
  )

  expect_type(out, "list")
  expect_true(all(c("subsetted", "sampled", "n_filtered", "filter_summary") %in% names(out)))

  expect_true(is.matrix(out$subsetted))
  expect_true(is.matrix(out$sampled))
  expect_true(is.numeric(out$n_filtered))
  expect_true(is.data.frame(out$filter_summary))

  # basic dimensions
  expect_equal(nrow(out$subsetted), length(inp$series.obs))
  expect_equal(nrow(out$sampled), length(inp$series.obs))
  expect_lte(ncol(out$sampled), 5)
})

test_that("filter_warm_simulations is deterministic given seed", {
  inp <- make_inputs(n_years = 30, n_real = 12, n_period = 10)

  out1 <- filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    power.obs = inp$power.obs,
    power.sim = inp$power.sim,
    power.period = inp$power.period,
    power.signif = inp$power.signif,
    sample.num = 5,
    seed = 999,
    save.plots = FALSE
  )

  out2 <- filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    power.obs = inp$power.obs,
    power.sim = inp$power.sim,
    power.period = inp$power.period,
    power.signif = inp$power.signif,
    sample.num = 5,
    seed = 999,
    save.plots = FALSE
  )

  expect_equal(out1$sampled, out2$sampled)
  expect_equal(out1$n_filtered, out2$n_filtered)
})

# -----------------------------------------------------------------------------
# 2) Input validation: isolate each fail-fast check
# -----------------------------------------------------------------------------

test_that("filter_warm_simulations validates series.obs vs series.sim year length", {
  inp <- make_inputs(n_years = 30, n_real = 10, n_period = 12)

  # break only year length: series.obs length != nrow(series.sim)
  series.obs_bad <- inp$series.obs[-1]  # length 29

  expect_error(
    filter_warm_simulations(
      series.obs = series.obs_bad,
      series.sim = inp$series.sim,
      power.obs = inp$power.obs,
      power.sim = inp$power.sim,
      power.period = inp$power.period,
      power.signif = inp$power.signif,
      save.plots = FALSE
    ),
    "Length mismatch: length\\(series\\.obs\\)"
  )
})

test_that("filter_warm_simulations validates power.period vs power.signif length", {
  inp <- make_inputs(n_years = 30, n_real = 10, n_period = 12)

  power.period_bad <- inp$power.period[-1]  # length 11, power.signif length 12

  expect_error(
    filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = inp$series.sim,
      power.obs = inp$power.obs,
      power.sim = inp$power.sim,
      power.period = power.period_bad,
      power.signif = inp$power.signif,
      save.plots = FALSE
    ),
    "must have the same length"
  )
})

test_that("filter_warm_simulations validates realization mismatch (ncol series.sim vs ncol power.sim)", {
  inp <- make_inputs(n_years = 30, n_real = 10, n_period = 12)

  # break only ncol(power.sim): keep all other dims consistent
  power.sim_bad <- inp$power.sim[, 1:2, drop = FALSE]  # now 2 realizations vs 10

  expect_error(
    filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = inp$series.sim,
      power.obs = inp$power.obs,
      power.sim = power.sim_bad,
      power.period = inp$power.period,
      power.signif = inp$power.signif,
      save.plots = FALSE
    ),
    "Realization mismatch: ncol\\(series\\.sim\\)"
  )
})

test_that("filter_warm_simulations validates period/grid mismatch (nrow power.sim)", {
  inp <- make_inputs(n_years = 30, n_real = 10, n_period = 12)

  # break only nrow(power.sim): keep ncol(power.sim) consistent
  power.sim_bad <- inp$power.sim[1:10, , drop = FALSE]  # 10 periods vs 12

  expect_error(
    filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = inp$series.sim,
      power.obs = inp$power.obs,
      power.sim = power.sim_bad,
      power.period = inp$power.period,
      power.signif = inp$power.signif,
      save.plots = FALSE
    ),
    "Period/grid mismatch: nrow\\(power\\.sim\\)"
  )
})

# -----------------------------------------------------------------------------
# 3) Priority 1.1 behavior: "most periods" power filter + fallback chain warnings
# -----------------------------------------------------------------------------

test_that("fallback activates when no realizations pass power filter (captures multiple warnings)", {
  inp <- make_inputs(n_years = 30, n_real = 10, n_period = 12)

  # Construct a situation where power filter is impossible to satisfy:
  # - Make observed significant at many periods
  # - Make simulated power always below significance so has_signal becomes FALSE for all realizations
  power.signif <- rep(10, length(inp$power.signif))
  power.obs <- rep(20, length(inp$power.obs))     # obs > signif => significant periods exist
  power.sim <- matrix(0.0, nrow = length(power.signif), ncol = ncol(inp$series.sim))  # never exceeds signif

  cap <- capture_warnings(
    filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = inp$series.sim,
      power.obs = power.obs,
      power.sim = power.sim,
      power.period = inp$power.period,
      power.signif = power.signif,
      sample.num = 5,
      seed = 123,
      save.plots = FALSE,
      bounds = list(
        mean = 0.0001, sd = 0.0001, min = 0.0001, max = 0.0001,
        sig.thr = 0.8, nsig.thr = 1.5,
        sig.frac = 0.8, nsig.frac = 0.95
      )
    )
  )

  out <- cap$value
  warns <- cap$warnings

  # Expect at least the initial fallback warning; subsequent warnings may occur depending on stats tightness.
  expect_true(any(grepl("Dropping power constraint", warns)))

  # Ensure function still returns a valid result
  expect_true(is.matrix(out$sampled))
  expect_equal(nrow(out$sampled), length(inp$series.obs))
  expect_lte(ncol(out$sampled), 5)

  # Ensure relaxation_level is reported (and not "none" in this constructed case)
  expect_true("relaxation_level" %in% names(out$filter_summary))
  expect_true(any(out$filter_summary$relaxation_level != "none"))
})

test_that("Priority 1.1 power filter respects bounds$sig.frac (80% default)", {
  inp <- make_inputs(n_years = 30, n_real = 10, n_period = 10)

  # Make 10 significant periods in obs by setting obs > signif everywhere
  power.signif <- rep(1, 10)
  power.obs <- rep(2, 10)

  # Construct power.sim so that realization 1 matches 8/10 periods (80%), realization 2 matches 7/10 (70%)
  power.sim <- matrix(0.5, nrow = 10, ncol = 10)

  # Ensure has_signal TRUE for both (at least one period exceeds signif)
  power.sim[1, 1] <- 1.2
  power.sim[1, 2] <- 1.2

  # With sig.thr=0.8, lower bound = obs*0.8 = 1.6; set "within bounds" to 1.7
  power.sim[1:8, 1] <- 1.7  # 80% within
  power.sim[1:7, 2] <- 1.7  # 70% within

  out <- filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    power.obs = power.obs,
    power.sim = power.sim,
    power.period = inp$power.period,
    power.signif = power.signif,
    sample.num = 5,
    seed = 1,
    save.plots = FALSE,
    bounds = list(
      mean = 1, sd = 1, min = 1, max = 1,  # stats wide open
      sig.thr = 0.8,
      nsig.thr = 1000,                     # effectively disable nonsig constraint
      sig.frac = 0.8,
      nsig.frac = 1.0                      # must be in (0,1]
    )
  )

  # With stats wide open and power filter active, we should retain >=1 realization.
  expect_gte(out$n_filtered, 1)
})

# -----------------------------------------------------------------------------
# 4) save.plots requires output.path
# -----------------------------------------------------------------------------

test_that("save.plots=TRUE requires output.path", {
  inp <- make_inputs()

  expect_error(
    filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = inp$series.sim,
      power.obs = inp$power.obs,
      power.sim = inp$power.sim,
      power.period = inp$power.period,
      power.signif = inp$power.signif,
      save.plots = TRUE,
      output.path = NULL
    ),
    "output\\.path"
  )
})
