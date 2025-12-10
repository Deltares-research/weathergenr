
library(ggplot2)
library(devtools)
library(weathergenr)
library(dplyr)
library(tidyr)
library(microbenchmark)

fpath <- "C:/Users/taner/WS/weathergenr-testdata/"
rdata <- paste0(fpath, "genarateWeatherSeries_ntoum_m1_rlz3_05122025.Rdata")
load(rdata)

daily.sim <- rlz_sample
daily.obs <- obs_sample
save.plots <- TRUE
variables <- variables
variable.labels = variables
variable.units <- NULL
realization.num <- rlz_num
wet.quantile <- 0.3 #mc.wet.quantile
extreme.quantile <- mc.extreme.quantile

timez <- Sys.time()

  out <- evaluateWegen(
    output.path = output.path,
    daily.sim = rlz_sample,
    daily.obs = obs_sample,
    save.plots = save.plots,
    variables = variables,
    realization.num = realization.num,
    wet.quantile = wet.quantile,
    extreme.quantile = extreme.quantile)

Sys.time() - timez
