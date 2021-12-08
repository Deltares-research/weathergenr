
#devtools::install_github("tanerumit/gridwegen")

# Packages needed
library(gridwegen)

# Path to output files
out_path <- "C:/wegentest/ntoum/"
nc_file <-  "ntoum.nc"
nc_dimnames <- list(x = "lon", y = "lat", time = "time")
variables <- c("precip", "temp", "temp_min", "temp_max")
variable_labels = c("Precipitation", "Mean Temperature", "Minimum Temperature", "Maximum Temperature")
variable_units = c("mm/day", "DegC", "DegC", "DegC")

# Read-in gridded weather data from netcdf
nc_data_ini <- readNetcdf(
    nc.path = paste0(out_path, "data/"),
    nc.file = nc_file,
    nc.dimnames = nc_dimnames,
    nc.variables = variables,
    origin.date = as.Date("1981-01-01"),
    has.leap.days = TRUE)

nc_data <- nc_data_ini

# Climate change perturbations
precip_changes <- list()
precip_changes$increments <- 2
precip_changes$mean <- list()
precip_changes$var  <- list()

temp_changes <- list()
temp_changes$increments <- 2
temp_changes$mean <- list()

############################ Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
precip_changes$mean$min <- c(0.7, 0.7, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
precip_changes$mean$max <- c(1.3, 1.3, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
precip_changes$var$min  <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
precip_changes$var$max  <- c(1.3, 1.3, 1.5, 1.5, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
temp_changes$mean$min   <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
temp_changes$mean$max   <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)

simulateWeather(
  climate.data = nc_data$data,
  climate.grid = nc_data$coords,
  year.start = 1981,
  year.num = 39,
  month.start = 1,
  variable.names = variables,
  variable.labels = variable_labels,
  variable.units = variable_units,
  sim.year.start = 2020,
  sim.year.num = 40,
  realization.num = 3,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.size = 20000,
  save.warm.results = TRUE,
  knn.annual.sample.size = 100,
  evaluate.model = FALSE,
  evaluate.grid.num = 30,
  apply.delta.changes = TRUE,
  apply.step.changes = TRUE,
  delta.precip = precip_changes,
  delta.temp = temp_changes,
  save.scenario.matrix = TRUE,
  output.path = out_path,
  output.ncfile.template = nc_data,
  output.ncfile.prefix = "clim_change_rlz"
)

