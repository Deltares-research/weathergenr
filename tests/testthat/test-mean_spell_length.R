

var <- c(
  9.12, 15.59, 4.84, 5.2, 1.54, 8.35, 14.37, 11.45,
  9.15, 0.81, 0.73, 1.05, 11.7, 12.35, 1.55, 3.9,
  7.19, 13.22, 8.87, 2.2, 9.24, 12.38, 14.69, 33.85,
  6.07, 0.68, 10.75, 9.25, 8.86, 3.4
)

testthat::test_that("Check avg dry spells (precip < 2 mm)", {
  testthat::expect_equal(mean_spell_length(x = var, threshold = 2, below = TRUE), 1.5)
})

testthat::test_that("Check avg wet spells (precip > 4 mm)", {
  testthat::expect_equal(mean_spell_length(x = var, threshold = 2, below = FALSE), 4.8)
})
