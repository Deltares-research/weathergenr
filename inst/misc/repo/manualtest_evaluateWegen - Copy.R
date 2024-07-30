

library(dplyr)
library(ggplot2)
library(weathergenr)
library(tidyr)

output_path <- "./inst/extdata/"
ncfile <- paste0(output_path, "ntoum_era5_data.nc")
ncdata <- readNetcdf(ncfile)

output <- generateWeatherSeries(
  weather.data = ncdata$data,
  weather.grid = ncdata$grid,
  weather.date = ncdata$date,
  variable.names = c("precip", "temp", "temp_min", "temp_max"),
  variable.labels = c("precip", "temp", "temp_min", "temp_max"),
  variable.units = NULL,
  sim.year.num = 30,
  sim.year.start = 2000,
  month.start = 1,
  realization.num = 5,
  warm.variable = "precip",
  warm.signif.level = 0.80,
  warm.sample.num = 10000,
  warm.subset.criteria = list(
    mean = c(0.90,1.10),
    sd = c(0.90,1.10),
    min = c(0.90,1.10),
    max = c(0.90,1.10),
    power = c(0.90, 10.00),
    nonsignif.threshold = 1.20),
  knn.sample.num = 100,
  mc.wet.quantile= 0.2,
  mc.extreme.quantile = 0.8,
  output.path = output_path,
  seed = 1111)

