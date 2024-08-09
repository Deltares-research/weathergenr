
weather.data = ncdata$data
weather.grid = ncdata$grid
weather.date = ncdata$date
variable.names = variables
variable.labels = variables
variable.units = NULL
sim.year.num = 52.0
sim.year.start = 2020
month.start = 1
realization.num = realization_num
warm.variable = "precip"
warm.signif.level = 0.80
warm.sample.num = 10000
warm.subset.criteria = NULL
knn.sample.num = 120
mc.wet.quantile= 0.2
mc.extreme.quantile = 0.8
evaluate.model = FALSE
evaluate.grid.num = 20
output.path = output_path
compute.parallel = FALSE
seed = 1


series.obs = warm_variable
series.sim = sim_annual
power.obs = warm_power$GWS
power.sim = sim_power
power.period = warm_power$GWS_period
power.signif = warm_power$GWS_signif
sample.num = realization.num
output.path = plots_path
bounds = warm.subset.criteria
seed = seed
save.series = FALSE
