
# Packages needed
library(gridwegen)

# Path to output files
path0 <- "C:/wegentest/"
out_path <- paste0(path0, "TEST/")

# Read-in gridded weather data from netcdf
nc_path <- system.file('extdata', package = 'gridwegen')
nc_file <- "ntoum.nc"
nc_dimnames <- list(x = "lon", y = "lat", time = "time")
nc_variables <- c("precip", "temp", "temp_min", "temp_max")
origin_date <- as.Date("1981-01-01")

nc_data <- readNetcdf(
    nc.path = nc_path,
    nc.file = nc_file,
    nc.dimnames = nc_dimnames,
    nc.variables = nc_variables,
    origin.date = origin_date,
    has.leap.days = TRUE)

# Simulate daily realizations of historical weather conditions
out <- simulateWeather(
  output.path = out_path,
  hist.climate = nc_data$data,
  hist.date.start = origin_date,
  variables = names(nc_data$variables)[1:4],
  sim.year.num = 40,
  sim.year.start = 2020,
  rlz.num = 5,
  return.date.indices = TRUE)



for(n in 1:rlz_num) {

  imposeClimateChanges(
      input.data = out[[n]],
      coordGrid = nc_data$tidy_data %>% dplyr::select(-data),
      file.suffix = n,
      output.path = paste0(out_path,"future/"),
      sim.date.start = sim_origin_date,
      variables = wg_variables,
      variable.units = wg_variable_units,
      nc.dimnames = nc_dimnames,
      change.settings = paste0(nc_path, "change_factors.xlsx"),
      save.scenario.matrix = FALSE,
      step = TRUE,
  )

}

