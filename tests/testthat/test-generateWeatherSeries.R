

testthat::test_that("check if the output dimensions are correct", {

  is.date <- function(x) inherits(x, 'Date')

  # Read weather data
  ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  ncdata <- readNetcdf(ncfile)
  rlz_num <- 3

  res <- generateWeatherSeries(
    weather.data = ncdata$data,
    weather.grid = ncdata$grid,
    weather.date = ncdata$date,
    variable.names = c("precip", "temp", "temp_min", "temp_max"),
    variable.labels = c("precip", "temp", "temp_min", "temp_max"),
    variable.units = NULL,
    sim.year.num = 20,
    sim.year.start = 2020,
    month.start = 1,
    realization.num = rlz_num,
    warm.variable = "precip",
    warm.signif.level = 0.90,
    warm.sample.num = 5000,
    warm.subset.criteria = list(mean = 0.1, sd = 0.2, min = 0.3,
       max = 0.3, sig.thr = 0.5, nsig.thr = 1.5),
    knn.sample.num = 100,
    mc.wet.quantile= 0.3,
    mc.extreme.quantile = 0.8,
    dry.spell.change = rep(1,12),
    wet.spell.change = rep(1,12),
    output.path = tempdir(),
    seed = 2024,
    compute.parallel = FALSE,
    num.cores = NULL,
    save.rdata = FALSE)

  res2 <- generateWeatherSeries(
    weather.data = ncdata$data,
    weather.grid = ncdata$grid,
    weather.date = ncdata$date,
    variable.names = c("precip", "temp", "temp_min", "temp_max"),
    variable.labels = c("precip", "temp", "temp_min", "temp_max"),
    variable.units = NULL,
    sim.year.num = 20,
    sim.year.start = 2020,
    month.start = 1,
    realization.num = rlz_num,
    warm.variable = "precip",
    warm.signif.level = 0.90,
    warm.sample.num = 5000,
    warm.subset.criteria = list(mean = 0.1, sd = 0.2, min = 0.3,
                                max = 0.3, sig.thr = 0.5, nsig.thr = 1.5),
    knn.sample.num = 100,
    mc.wet.quantile= 0.3,
    mc.extreme.quantile = 0.8,
    dry.spell.change = rep(1,12),
    wet.spell.change = rep(1,12),
    output.path = tempdir(),
    seed = 2024,
    compute.parallel = FALSE,
    num.cores = NULL,
    save.rdata = FALSE)

  res3 <- generateWeatherSeries(
    weather.data = ncdata$data,
    weather.grid = ncdata$grid,
    weather.date = ncdata$date,
    variable.names = c("precip", "temp", "temp_min", "temp_max"),
    variable.labels = c("precip", "temp", "temp_min", "temp_max"),
    variable.units = NULL,
    sim.year.num = 20,
    sim.year.start = 2020,
    month.start = 1,
    realization.num = rlz_num,
    warm.variable = "precip",
    warm.signif.level = 0.90,
    warm.sample.num = 5000,
    warm.subset.criteria = list(mean = 0.1, sd = 0.2, min = 0.3,
                                max = 0.3, sig.thr = 0.5, nsig.thr = 1.5),
    knn.sample.num = 100,
    mc.wet.quantile= 0.3,
    mc.extreme.quantile = 0.8,
    dry.spell.change = rep(1,12),
    wet.spell.change = rep(1,12),
    output.path = tempdir(),
    seed = 9999,
    compute.parallel = FALSE,
    num.cores = NULL,
    save.rdata = FALSE)

  # Structure checks
  expect_type(res, "list")
  expect_true(all(c("resampled", "dates") %in% names(res)))
  expect_s3_class(res$resampled, "tbl_df")
  expect_true(is.date(res$dates))
  expect_equal(ncol(res$resampled), rlz_num)  # matches realization.num
  expect_equal(length(res$dates), nrow(res$resampled))

  # Check stochastic behavior
  expect_equal(res, res2)
  expect_failure(expect_equal(res, res3))

})
