

testthat::test_that("check if the output dimensions are correct", {

  # Read weather data
  ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  ncdata <- readNetcdf(ncfile)

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
    realization.num = 5,
    warm.variable = "precip",
    warm.signif.level = 0.90,
    warm.sample.num = 5000,
    warm.subset.criteria = list(mean = 0.1, sd = 0.2, min = 0.3,
       max = 0.3, signif.threshold = 0.5, nonsignif.threshold = 1.5),
    knn.sample.num = 100,
    mc.wet.quantile= 0.3,
    mc.extreme.quantile = 0.8,
    output.path = tempdir(),
    seed = 2024,
    compute.parallel = FALSE,
    num.cores = NULL,
    save.rdata = FALSE)

  resampled_dates_ini <- read.csv("../testdata/resampled_dates.csv")
  resampled_dates <- dplyr::bind_cols(lapply(resampled_dates_ini, as.Date))

  expect_equal(resampled_dates, output$resampled, ignore_attr = TRUE)

})
