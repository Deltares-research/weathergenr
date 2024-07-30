

library(dplyr)
library(ggplot2)
library(weathergenr)
library(tidyr)

output_path <- "C:/Users/taner/Workspace/test-weathergenr/ntoum_test/"

ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
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
  realization.num = 11,
  warm.variable = "precip",
  warm.signif.level = 0.80,
  warm.sample.num = 1000,
  warm.subset.criteria = list(
    mean = c(0.90,1.10),
    sd = c(0.80,1.20),
    min = c(0.80,1.20),
    max = c(0.80,1.20),
    power = c(0.20, 10.00),
    nonsignif.threshold = 1.20),
  knn.sample.num = 100,
  mc.wet.quantile= 0.2,
  mc.extreme.quantile = 0.8,
  output.path = output_path,
  seed = 1111)

#########################################

#Comparisons
resampled_dates <- readr::read_csv(paste0(output_path, "resampled_dates.csv"))
simulated_dates <- readr::read_csv(paste0(output_path, "sim_dates.csv"))
realization.num <- ncol(resampled_dates)

evaluate.grid.num <- 20

grids <- ncdata$grid$id
sampleGrids <- sample(grids, size = min(evaluate.grid.num, length(grids)))

daily_obs <- lapply(ncdata$data[sampleGrids], function(x)
  dplyr::mutate(x, date = ncdata$date, .before = 1))

daily_sim <- vector(mode = "list", length = realization.num)

for (n in 1:realization.num) {

  # Resample order
  day_order <- match(resampled_dates[[n]], ncdata$date)

  # Obtain stochastic series by re-ordering historical data
  daily_sim[[n]] <- lapply(ncdata$data[sampleGrids], function(x)
    x[day_order,] %>% mutate(date = simulated_dates[[1]], .before = 1))
}


evaluateWegen(daily.sim = daily_sim,
              daily.obs = daily_obs,
              output.path = paste0(output_path, "plots"),
              variables = c("precip", "temp", "temp_min", "temp_max"),
              variable.labels = c("precip", "temp", "temp_min", "temp_max"),
              variable.units = NULL,
              realization.num = realization.num,
              wet.quantile = 0.2,
              extreme.quantile = 0.8,
              show.title = TRUE)














load("../test-weathergenr/input_evaluate_wegen.Rdata")

wet.quantile = 0.3
extreme.quantile = 0.8
daily.sim = rlz_sample
daily.obs = obs_sample
output.path = plots_path
variables = variable.names
variable.labels = variable.names
variable.units = NULL
realization.num = realization.num
show.title = TRUE

evaluateWegen(daily.sim = rlz_sample,
              daily.obs = obs_sample,
              output.path = paste0(plots_path,"plots"),
              variables = variable.names,
              variable.labels = variable.names,
              variable.units = variable.units,
              realization.num = realization.num,
              wet.quantile = wet.quantile,
              extreme.quantile = extreme.quantile,
              show.title = TRUE)


