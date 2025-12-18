library(testthat)

# ------------------------------------------------------------------------------
test_that("filter_warm_simulations validates inputs", {

  series.obs <- rnorm(10)
  series.sim <- matrix(rnorm(100), nrow = 10)
  power.obs <- runif(5)
  power.sim <- matrix(runif(50), nrow = 5)
  power.period <- 1:5
  power.signif <- rep(0.5, 5)

  expect_error(
    filter_warm_simulations(series.obs = "a"),
    "'series.obs' must be a numeric vector"
  )

  expect_error(
    filter_warm_simulations(series.obs, series.sim = 1:10),
    "'series.sim' must be a numeric matrix"
  )

  expect_error(
    filter_warm_simulations(series.obs, series.sim, power.obs = "x"),
    "'power.obs' must be a numeric vector"
  )

  expect_error(
    filter_warm_simulations(series.obs, series.sim,
                            power.obs, power.sim,
                            power.period = 1:3,
                            power.signif = 1:5),
    "same length"
  )

  expect_error(
    filter_warm_simulations(series.obs, series.sim,
                            power.obs, power.sim[, 1:2],
                            power.period, power.signif),
    "same number of realizations"
  )
})

# ------------------------------------------------------------------------------
test_that("filter_warm_simulations runs and returns correct structure", {

  set.seed(1)

  series.obs <- rnorm(20)
  series.sim <- replicate(10, series.obs + rnorm(20, sd = 0.1))

  power.obs <- runif(6, 1, 2)
  power.sim <- matrix(replicate(10, power.obs), nrow = 6)

  res <- suppressWarnings(
    filter_warm_simulations(
      series.obs,
      series.sim,
      power.obs,
      power.sim,
      power.period = 1:6,
      power.signif = rep(0.8, 6),
      save.plots = FALSE,
      save.series = FALSE,
      seed = 42
    )
  )

  expect_named(res, c("subsetted", "sampled", "n_filtered", "filter_summary"))
  expect_true(is.matrix(res$subsetted))
  expect_true(is.matrix(res$sampled))
  expect_true(is.numeric(res$n_filtered))
  expect_true(is.data.frame(res$filter_summary))
})

# ------------------------------------------------------------------------------
test_that("statistical filters remove bad realizations", {

  set.seed(2)

  series.obs <- rnorm(30)
  good <- replicate(5, series.obs + rnorm(30, sd = 0.05))
  bad  <- replicate(5, series.obs * 5)
  series.sim <- cbind(good, bad)

  power.obs <- runif(4, 1, 2)
  power.sim <- matrix(rep(power.obs, 10), nrow = 4)

  res <- suppressWarnings(
    filter_warm_simulations(
      series.obs,
      series.sim,
      power.obs,
      power.sim,
      power.period = 1:4,
      power.signif = rep(0.5, 4),
      save.plots = FALSE,
      save.series = FALSE
    )
  )

  expect_lt(res$n_filtered, ncol(series.sim))
  expect_gt(res$n_filtered, 0)
})

# ------------------------------------------------------------------------------
test_that("power spectrum filter excludes realizations without signal", {

  series.obs <- rnorm(20)
  series.sim <- replicate(6, series.obs + rnorm(20, 0.1))

  power.obs <- c(2, 1, 0.5)
  power.signif <- c(1, 1, 1)

  power.sim <- cbind(
    matrix(c(2, 1.2, 0.6), nrow = 3, ncol = 3),
    matrix(c(0.3, 0.2, 0.1), nrow = 3, ncol = 3)
  )

  res <- suppressWarnings(
    filter_warm_simulations(
      series.obs,
      series.sim,
      power.obs,
      power.sim,
      power.period = 1:3,
      power.signif = power.signif,
      save.plots = FALSE,
      save.series = FALSE
    )
  )

  expect_gt(res$n_filtered, 0)
  expect_lte(res$n_filtered, ncol(series.sim))
})

# ------------------------------------------------------------------------------
test_that("fallback activates when no realization passes power filter", {

  series.obs <- rnorm(15)
  series.sim <- replicate(5, series.obs + rnorm(15, 0.1))

  power.obs <- c(10, 10)
  power.sim <- matrix(0.01, nrow = 2, ncol = 5)

  expect_warning(
    res <- filter_warm_simulations(
      series.obs,
      series.sim,
      power.obs,
      power.sim,
      power.period = 1:2,
      power.signif = c(5, 5),
      save.plots = FALSE,
      save.series = FALSE
    ),
    regexp = ".*"
  )

  expect_gt(res$n_filtered, 0)
})

# ------------------------------------------------------------------------------
test_that("sampling is reproducible with seed", {

  series.obs <- rnorm(20)
  series.sim <- replicate(20, series.obs + rnorm(20, 0.1))

  power.obs <- runif(5, 1, 2)
  power.sim <- matrix(rep(power.obs, 20), nrow = 5)

  res1 <- suppressWarnings(
    filter_warm_simulations(
      series.obs, series.sim,
      power.obs, power.sim,
      power.period = 1:5,
      power.signif = rep(0.5, 5),
      seed = 123,
      save.plots = FALSE,
      save.series = FALSE
    )
  )

  res2 <- suppressWarnings(
    filter_warm_simulations(
      series.obs, series.sim,
      power.obs, power.sim,
      power.period = 1:5,
      power.signif = rep(0.5, 5),
      seed = 123,
      save.plots = FALSE,
      save.series = FALSE
    )
  )

  expect_equal(res1$sampled, res2$sampled)
})

# ------------------------------------------------------------------------------
test_that("filter_summary is consistent", {

  series.obs <- rnorm(10)
  series.sim <- replicate(5, series.obs)

  power.obs <- runif(3, 1, 2)
  power.sim <- matrix(rep(power.obs, 5), nrow = 3)

  res <- suppressWarnings(
    filter_warm_simulations(
      series.obs,
      series.sim,
      power.obs,
      power.sim,
      power.period = 1:3,
      power.signif = rep(0.5, 3),
      save.plots = FALSE,
      save.series = FALSE
    )
  )

  expect_equal(
    res$filter_summary$n_passed[res$filter_summary$filter == "mean"],
    5
  )
})
