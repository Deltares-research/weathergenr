

library(dplyr)
library(ggplot2)
library(weathergenr)
library(tidyr)

load("./tests/input_evaluate_wegen.Rdata")
source("./R/evaluateWegen.R")


evaluateWegen(daily.sim = rlz_sample,
              daily.obs = obs_sample,
              output.path = plots_path,
              variables = variable.names,
              variable.labels = variable.labels,
              variable.units = variable.units,
              realization.num = realization.num,
              wet.quantile = mc.wet.quantile,
              extreme.quantile = mc.extreme.quantile,
              show.title = TRUE)


daily.sim = rlz_sample
daily.obs = obs_sample
output.path = plots_path
variables = variable.names
variable.labels = variable.labels
variable.units = variable.units
realization.num = realization.num
wet.quantile = mc.wet.quantile
extreme.quantile = mc.extreme.quantile
show.title = TRUE
