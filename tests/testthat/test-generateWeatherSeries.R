

testthat::test_that("check if the output dimensions are correct", {

 sim_year_start = 2020
 sim_year_num = 20
 month_start = 1
 realization_num = 2

 date_begin <- as.Date(paste(sim_year_start, month_start, "01", sep = "-"))
 date_end <- as.Date(paste(sim_year_start+sim_year_num, month_start, "01", sep = "-"))-1
 date_vector <- seq.Date(date_begin, date_end, by = "day")

 output <- generateWeatherSeries(
    weather.data = ncdata$data,
    weather.grid = ncdata$grid,
    weather.date = ncdata$date,
    variable.names = c("precip", "temp", "temp_min", "temp_max"),
    variable.labels = c("precip", "temp", "temp_min", "temp_max"),
    variable.units = NULL,
    sim.year.num = sim_year_num,
    sim.year.start = sim_year_start,
    month.start = month_start,
    realization.num = realization_num,
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

  # Check if the resampled dates are from the historic days
  expect_true(all(sapply(1:realization_num, function(x) output$resampled[[x]] %in% ncdata$date)))

  # Check future date vector
  expect_equal(output$dates, date_vect)

})
