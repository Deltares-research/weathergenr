
library(devtools)
library(weathergenr)
library(dplyr)
library(tidyr)

output_path <- "C:/TEMP/TEST8"
variables <- c("precip", "temp", "temp_min", "temp_max")
variable_labels <- c("Precip.", "Temp. (avg)", "Temp. (min)", "Temp. (max)")
realization_num <- 5

ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
ncdata <- readNetcdf(ncfile)

mc.wet.quantile= 0.2
mc.extreme.quantile = 0.8
warm.subset.criteria = list(mean = 0.05, sd = 0.05, min = 0.05,
                            max = 0.05, sig.thr = 0.8, nsig.thr = 1.5)

weather.data = ncdata$data
weather.grid = ncdata$grid
weather.date = ncdata$date
variable.names = variables
variable.labels = variables
sim.year.num = 20
sim.year.start = 2020
month.start = 1
realization.num = realization_num
warm.variable = "precip"
warm.signif.level = 0.90
warm.sample.num = 20000
warm.subset.criteria = warm.subset.criteria
knn.sample.num = 100
mc.wet.quantile= mc.wet.quantile
mc.extreme.quantile = mc.extreme.quantile
dry.spell.change = rep(1, 12)
wet.spell.change = rep(1, 12)
output.path = output_path
compute.parallel = FALSE
num.cores = NULL
save.rdata = FALSE
seed = 143434

stochastic_weather <- generateWeatherSeries(
  weather.data = ncdata$data,
  weather.grid = ncdata$grid,
  weather.date = ncdata$date,
  variable.names = variables,
  variable.labels = variables,
  sim.year.num = 20,
  sim.year.start = 2020,
  month.start = 1,
  realization.num = realization_num,
  warm.variable = "precip",
  warm.signif.level = 0.70,
  warm.sample.num = 20000,
  warm.subset.criteria = warm.subset.criteria,
  knn.sample.num = 100,
  mc.wet.quantile= mc.wet.quantile,
  mc.extreme.quantile = mc.extreme.quantile,
  dry.spell.change = rep(1, 12),
  wet.spell.change = rep(1, 12),
  output.path = output_path,
  compute.parallel = FALSE,
  num.cores = NULL,
  seed = 143434)

day_order <- sapply(1:realization_num,
                    function(n) match(stochastic_weather$resampled[[n]], ncdata$date))

rlz_sample <- list()
for (n in 1:realization_num) {
  rlz_sample[[n]] <- lapply(ncdata$data[ncdata$grid$id], function(x) x[day_order[,n],] %>%
                              select(precip,  temp, temp_min, temp_max) %>%
                              mutate(date = stochastic_weather$dates, .before = 1))
}

obs_sample <- lapply(ncdata$data[ncdata$grid$id], function(x) x %>%
                       select(precip, temp, temp_min, temp_max) %>%
                       dplyr::mutate(date = ncdata$date, .before = 1))

out <- evaluateWegen(daily.sim = rlz_sample,
                     daily.obs = obs_sample,
                     output.path = output_path,
                     save.plots = TRUE,
                     variables = variables,
                     variable.labels = variable_labels,
                     variable.units = NULL,
                     realization.num = realization_num,
                     wet.quantile = mc.wet.quantile,
                     extreme.quantile = mc.extreme.quantile)
