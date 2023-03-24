


library(weathergenr)
output_path <- tempdir()
variables <- c("precip", "temp", "temp_min", "temp_max")
realization_num <- 3

## ----stochastic2, results='hide', eval = TRUE, cache=TRUE---------------------
stochastic_weather <- generateWeatherSeries(
  weather.data = ncdata$data,
  weather.grid = ncdata$grid,
  weather.date = ncdata$date,
  variable.names = variables,
  variable.labels = variables,
  variable.units = NULL,
  sim.year.num = 20,
  sim.year.start = 2020,
  month.start = 1,
  realization.num = realization_num,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.num = 10000,
  warm.subset.criteria = NULL,
  knn.sample.num = 120,
  mc.wet.quantile= 0.2,
  mc.extreme.quantile = 0.8,
  evaluate.model = FALSE,
  evaluate.grid.num = 20,
  output.path = output_path,
  seed = 123)

