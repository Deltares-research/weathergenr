



load("../test-weathergenr/testdata_waveletARsubset.rdata")
library(dplyr)
library(ggplot2)
library(tidyr)

warm_criteria = list(
  mean = c(0.90,1.10),
  sd = c(0.80,1.20),
  min = c(0.80,1.20),
  max = c(0.80,1.20),
  power = c(0.20,10),
  nonsignif.threshold = 1.20)

# subsetting from generated warm series
sim_annual_sub <- waveletARSubset(
  series.obs = warm_variable,
  series.sim = sim_annual,
  power.obs = warm_power$GWS,
  power.sim = sim_power,
  power.period = warm_power$GWS_period,
  power.signif = warm_power$GWS_signif,
  sample.num = realization.num,
  output.path = plots_path,
  bounds = warm_criteria,
  seed = seed,
  save.series = FALSE)


series.obs = warm_variable
series.sim = sim_annual
power.obs = warm_power$GWS
power.sim = sim_power
power.period = warm_power$GWS_period
power.signif = warm_power$GWS_signif
sample.num = realization.num
seed = seed
save.plots = TRUE
save.series = FALSE
output.path = plots_path
padding = TRUE
bounds = warm_criteria
