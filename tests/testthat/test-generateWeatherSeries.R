#
# make_long_weather_data <- function(start = "1980-01-01", end = "2004-12-31") {
#
#   dates <- seq.Date(as.Date(start), as.Date(end), by = "day")
#   dates <- dates[format(dates, "%m-%d") != "02-29"]
#
#   n <- length(dates)
#
#   list(
#     dates = dates,
#     weather.data = list(
#       data.frame(
#         precip = rgamma(n, shape = 2, scale = 2),
#         temp   = rnorm(n, mean = 10, sd = 5)
#       )
#     ),
#     weather.grid = data.frame(
#       id = 1,
#       xind = 1,
#       yind = 1,
#       x = 0,
#       y = 0
#     )
#   )
# }
#
#
# test_that("generateWeatherSeries returns correct structure with long data", {
#
#   set.seed(123)
#
#   d <- make_long_weather_data()
#
#   out <- generateWeatherSeries(
#     weather.data = d$weather.data,
#     weather.grid = d$weather.grid,
#     weather.date = d$dates,
#     variables = c("precip", "temp"),
#     sim.year.num = 5,
#     sim.year.start = 2020,
#     realization.num = 3,
#     month.start = 1,
#     warm.sample.num = 200,
#     knn.sample.num = 20,
#     seed = 42,
#     compute.parallel = FALSE,
#     output.path = tempdir()
#   )
#
#   expect_type(out, "list")
#   expect_true(all(c("resampled", "dates") %in% names(out)))
#
#   expect_s3_class(out$dates, "Date")
#   expect_equal(nrow(out$resampled), length(out$dates))
#   expect_equal(ncol(out$resampled), 3)
#   expect_false(any(is.na(out$dates)))
# })
#
# test_that("calendar-year simulation starts Jan 1 and has 365-day years", {
#
#   set.seed(1)
#
#   d <- make_long_weather_data()
#
#   out <- generateWeatherSeries(
#     weather.data = d$weather.data,
#     weather.grid = d$weather.grid,
#     weather.date = d$dates,
#     variables = c("precip", "temp"),
#     sim.year.num = 3,
#     sim.year.start = 2015,
#     realization.num = 1,
#     month.start = 1,
#     warm.sample.num = 200,
#     knn.sample.num = 20,
#     seed = 99,
#     output.path = tempdir()
#   )
#
#   dates <- out$dates
#
#   expect_equal(format(dates[1], "%m-%d"), "01-01")
#   expect_equal(length(dates), 3 * 365)
#
#   yrs <- split(dates, format(dates, "%Y"))
#   expect_true(all(lengths(yrs) == 365))
# })
#
#
# test_that("calendar-year simulation starts Jan 1 and has 365-day years", {
#
#   set.seed(1)
#
#   d <- make_long_weather_data()
#
#   out <- generateWeatherSeries(
#     weather.data = d$weather.data,
#     weather.grid = d$weather.grid,
#     weather.date = d$dates,
#     variables = c("precip", "temp"),
#     sim.year.num = 3,
#     sim.year.start = 2015,
#     realization.num = 1,
#     month.start = 1,
#     warm.sample.num = 200,
#     knn.sample.num = 20,
#     seed = 99,
#     output.path = tempdir()
#   )
#
#   dates <- out$dates
#
#   expect_equal(format(dates[1], "%m-%d"), "01-01")
#   expect_equal(length(dates), 3 * 365)
#
#   yrs <- split(dates, format(dates, "%Y"))
#   expect_true(all(lengths(yrs) == 365))
# })
#
# test_that("water-year simulation starts at month.start and spans 365 days", {
#
#   set.seed(2)
#
#   d <- make_long_weather_data()
#
#   out <- generateWeatherSeries(
#     weather.data = d$weather.data,
#     weather.grid = d$weather.grid,
#     weather.date = d$dates,
#     variables = c("precip", "temp"),
#     sim.year.num = 3,
#     sim.year.start = 2015,
#     realization.num = 1,
#     month.start = 10,
#     warm.sample.num = 200,
#     knn.sample.num = 20,
#     seed = 11,
#     output.path = tempdir()
#   )
#
#   dates <- out$dates
#
#   expect_equal(format(dates[1], "%m-%d"), "10-01")
#   expect_equal(length(dates), 3 * 365)
#
#   wyears <- get_water_year(dates, month.start = 10)
#   expect_true(all(table(wyears) == 365))
# })
#
# test_that("simulated dates never include February 29", {
#
#   set.seed(3)
#
#   d <- make_long_weather_data()
#
#   out <- generateWeatherSeries(
#     weather.data = d$weather.data,
#     weather.grid = d$weather.grid,
#     weather.date = d$dates,
#     variables = c("precip", "temp"),
#     sim.year.num = 4,
#     sim.year.start = 2018,
#     realization.num = 1,
#     month.start = 1,
#     warm.sample.num = 200,
#     knn.sample.num = 20,
#     seed = 7,
#     output.path = tempdir()
#   )
#
#   expect_false(any(format(out$dates, "%m-%d") == "02-29"))
# })
