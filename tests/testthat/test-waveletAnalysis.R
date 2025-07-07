test_that("waveletAnalysis returns expected structure and types", {
  set.seed(42)
  n <- 64
  # Simulated time series: sine with noise
  t <- 1:n
  y <- sin(2 * pi * t / 8) + rnorm(n, 0, 0.5)

  result <- waveletAnalysis(y, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("GWS", "GWS_signif", "GWS_period", "signif_periods"))
  expect_true(is.numeric(result$GWS))
  expect_true(is.numeric(result$GWS_signif))
  expect_true(is.numeric(result$GWS_period))
  expect_true(is.numeric(result$signif_periods) || is.integer(result$signif_periods))
  expect_length(result$GWS, length(result$GWS_signif))
  expect_length(result$GWS, length(result$GWS_period))
})

test_that("Significant periods are detected for strong periodic signal", {
  n <- 64
  t <- 1:n
  # Pure periodic signal (no noise)
  y <- sin(2 * pi * t / 8)
  result <- waveletAnalysis(y, plot = FALSE)
  # The expected period is 8
  expect_true(any(abs(result$GWS_period[result$signif_periods] - 8) < 2))
})

test_that("Function errors on NA or non-numeric input", {
  y <- 1:20
  y_na <- y
  y_na[5] <- NA
  expect_error(waveletAnalysis(y_na, plot = FALSE), "contains missing values")
  expect_error(waveletAnalysis(letters, plot = FALSE))
})

test_that("Function errors for invalid noise.type or signif.level", {
  y <- rnorm(16)
  expect_error(waveletAnalysis(y, noise.type = "invalid"), "noise.type")
  expect_error(waveletAnalysis(y, signif.level = -0.1), "signif.level")
  expect_error(waveletAnalysis(y, signif.level = 1.1), "signif.level")
})

test_that("Function detects no significant periodicity for white noise", {
  set.seed(99)
  y <- rnorm(64)
  result <- waveletAnalysis(y, signif.level = 0.95, plot = FALSE)
  # Should have zero or very few significant periods
  expect_true(length(result$signif_periods) <= 2)
})
