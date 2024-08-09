

testthat::test_that("check if the output dimensions are correct", {

  output <- generateWeatherSeries(
    weather.data = ncdata$data,
    weather.grid = ncdata$grid,
    weather.date = ncdata$date,
    variable.names = c("precip", "temp", "temp_min", "temp_max"),
    variable.labels = c("precip", "temp", "temp_min", "temp_max"),
    variable.units = NULL,
    sim.year.num = 20,
    sim.year.start = 2020,
    month.start = 1,
    realization.num = 2,
    warm.variable = "precip",
    warm.signif.level = 0.90,
    warm.sample.num = 1000,
    warm.subset.criteria = NULL,
    knn.sample.num = 100,
    mc.wet.quantile= 0.2,
    mc.extreme.quantile = 0.8,
    evaluate.model = FALSE,
    evaluate.grid.num = 20,
    output.path = tempdir(),
    seed = 123)

  expect_equal(nrow(output[[1]]), 20 * 365)
  expect_equal(length(output[[2]]), 20 * 365)

})
