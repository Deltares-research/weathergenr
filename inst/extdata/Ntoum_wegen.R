
# Todo
# Add dry/wet spell plots for performance
# Add step changes......

# Load package
#devtools::install_github("tanerumit/gridwegen@dev")

# Packages needed
library(ncdf4)
library(dplyr)
library(tidyr)
library(fitdistrplus)
library(gridwegen)

# Path to output files
path0 <- "C:/Users/taner/OneDrive - Stichting Deltares/_DELTARES/02 Projects/11206634 Gabon/05 Models/wegen/"
out_path <- paste0(path0, "TEST03/")
nc_path <- paste0(path0, "input/")

# Path to historical gridded data
nc_file <- "localP_chirpsP_era5T_19810101_chunks.nc"
nc_dimnames <- list(x = "lon", y = "lat", time = "time")
origin_date <- as.Date("1981-01-01")
sim_origin_date <- as.Date("2020-01-01")
wg_variables <- c("precip", "temp", "temp_min", "temp_max", "pet")
wg_variable_labels <- c("Precipitation", "Avg. Temperature", "Min. Temperature","Max. Temperature")
wg_variable_units  <- c("mm/day", "°C", "°C", "°C")

warm_variable <- "precip"

water_year_months <- 1:12
rlz_num = 5
sim_years = 40
ymax = 40
knn_annual_sample_size = 20
n=1

# Read-in gridded weather data from netcdf
nc_data <- readNetcdf(
    nc.path = nc_path,
    nc.file = nc_file,
    nc.dimnames = nc_dimnames,
    nc.variables = wg_variables,
    origin.date = origin_date,
    leap.year = TRUE)

out <- simulateWeather(
  output.path = out_path,
  hist.climate = nc_data$tidy_data,
  hist.date.start = origin_date,
  sim.year.start = as.numeric(format(sim_origin_date,"%Y")),
  month.list = water_year_months,
  variables = wg_variables,
  variable.labels = wg_variable_labels,
  variable.units = wg_variable_units,
  warm.variable = warm_variable,
  knn.annual.sample.size = knn_annual_sample_size,
  warm.signif.level = 0.90,
  ymax = sim_years,
  nmax = rlz_num)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


for(n in 1:nmax) {

  imposeClimateChanges(
      input.path = paste0(out_path,"historical/"),
      input.data = nvar_filenames[n],
      file.suffix = n,
      output.path = paste0(out_path,"future/"),
      sim.date.start = as.Date("2020-01-01"),
      variables = c("precip", "temp", "temp_min", "temp_max"),
      variable.units = c("mm/day", "°C", "°C", "°C"),
      nc.dimnames = list(x = "lon", y = "lat", time = "time"),
      change.settings = paste0(nc_path, "change_factors.xlsx"),
      save.scenario.matrix = FALSE,
      step = TRUE,
  )

}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
