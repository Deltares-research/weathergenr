
################################################################################
############## PERFORMANCE PROFILING ###########################################
################################################################################

# Libraries

library(profvis)
library(devtools)
library(weathergenr)
library(dplyr)
library(tidyr)
library(microbenchmark)


# Read-in forcing data & geometry
month.start <- 1

ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
ncdata <- read_netcdf(ncfile)

#### Define all variables in advance for testing
mc.wet.quantile <- 0.1
weather.data <- ncdata$data
weather.grid <- ncdata$grid
weather.date <- ncdata$date
variables <- c("precip", "temp", "temp_min", "temp_max")
variable.labels <- variables
sim.year.num <- 35
sim.year.start <- 2020
realization.num <- 5
warm.variable <- "precip"
warm.signif.level <- 0.90
warm.sample.num <- 10000
warm.subset.criteria <- list(mean = 0.05, sd = 0.05, min = 0.05, max = 0.05, sig.thr = 0.8, nsig.thr = 1.5)
knn.sample.num <- 100

mc.extreme.quantile <- 0.8
dry.spell.change <- rep(1, 12)
wet.spell.change <- rep(1, 12)
output.path <- output_path
seed <- 1242
compute.parallel<- FALSE
num.cores <- NULL

profvis({

  stochastic_weather <- generateWeatherSeries(
    weather.data = weather.data,
    weather.grid = weather.grid,
    weather.date = weather.date,
    variables = variables,
    sim.year.num = sim.year.num,
    sim.year.start = sim.year.start,
    month.start = month.start,
    realization.num = realization.num,
    warm.variable = warm.variable,
    warm.signif.level = warm.signif.level,
    warm.sample.num = warm.sample.num,
    warm.subset.criteria = warm.subset.criteria,
    knn.sample.num = knn.sample.num,
    mc.wet.quantile= mc.wet.quantile,
    mc.extreme.quantile = mc.extreme.quantile,
    dry.spell.change = dry.spell.change,
    wet.spell.change = wet.spell.change,
    output.path = output.path,
    compute.parallel = FALSE,
    num.cores = num.cores,
    seed = seed)


})
