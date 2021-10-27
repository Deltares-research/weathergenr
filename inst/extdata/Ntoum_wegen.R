

# Load package
#devtools::install_github("tanerumit/gridwegen@dev")

# Packages needed
library(gridwegen)
library(tidyr)
library(dplyr)

# Path to output files
path0 <- "C:/Users/taner/OneDrive - Stichting Deltares/_DELTARES/02 Projects/11206634 Gabon/05 Models/wegen/"

# Path to historical gridded data
out_path <- paste0(path0, "TEST02/")
nc_path <- paste0(path0, "input/")
nc_file <- "localP_chirpsP_era5T_19810101_chunks.nc"
nc_dimnames <- list(x = "lon", y = "lat", time = "time")
origin_date <- as.Date("1981-01-01")
sim_origin_date <- as.Date("2020-01-01")
wg_variables <- c("precip", "temp", "temp_min", "temp_max")
wg_variable_labels <- c("Precipitation", "Avg. Temperature", "Min. Temperature","Max. Temperature")
wg_variable_units  <- c("mm/day", "°C", "°C", "°C", "mm/day")
project_name = "ntoum_test"
month_list <- 1:12
rlz_num = 5
sim_years = 40


# Read-in gridded weather data from netcdf
nc_data <- readNetcdf(
    nc.path = nc_path,
    nc.file = nc_file,
    nc.dimnames = nc_dimnames,
    nc.variables = wg_variables,
    origin.date = origin_date,
    leap.year = TRUE)

sim_dates <- tibble(date = sim_origin_date + 0:(ymax*366)) %>%
  mutate(year = as.numeric(format(date,"%Y")),
         month = as.numeric(format(date,"%m")),
         day = as.numeric(format(date,"%d"))) %>%
  filter(month!=2 | day!=29) %>%
  slice(1:(ymax*365))

dates_res <- simulateWeather(
  output.path = out_path,
  hist.climate = nc_data$tidy_data,
  hist.date.start = origin_date,
  sim.year.start = as.numeric(format(sim_origin_date,"%Y")),
  month.list = month_list,
  variables = wg_variables,
  variable.labels = wg_variable_labels,
  variable.units = wg_variable_units,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  ymax = sim_years,
  nmax = rlz_num)

#### Calculate PET and save to netcdf file
for (n in 1:rlz_num) {

  day_order <- match(dates_res[[n]], nc_data$nc_dates)

  # Sample realizations
  rlz <- lapply(nc_data$tidy_data, function(x) x[day_order,])

  # Calculate PET based on Hargreaves Eq.
  rlz <- lapply(1:length(rlz), function(x)
    mutate(rlz[[x]], pet = hargreavesPet(months = sim_dates$month, temp = rlz[[x]]$temp,
            tdiff = rlz[[x]]$temp_max - rlz[[x]]$temp_min, lat = climate_tidy$y[x])))

  writeNetcdf(
      data = rlz,
      output.path = paste0(out_path,"historical/"),
      nc.dimensions = names(nc_data$nc_dimensions),
      nc.dimnames = nc_dimnames,
      origin.date = sim_origin_date,
      calendar.type = "no leap",
      variables = c(wg_variables, "pet"),
      variable.units = c(wg_variable_units, "mm/day"),
      file.suffix = n
  )
}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


# natural variability realizations in the input folder
nvar_filenames <- list.files(paste0(out_path,"historical/"))
nmax <- length(nvar_filenames)

for(n in 1:rlz_num) {

  imposeClimateChanges(
      proj.name = "ntoum",
      in.path = paste0(out_path,"historical/"),
      in.file = nvar_filenames[n],
      file.suffix = n,
      out.path = paste0(out_path,"future2/"),
      sim.date.begin = as.Date("2020-01-01"),
      wg.vars = c("precip", "temp", "temp_min", "temp_max"),
      wg.var.units = c("mm/day", "°C", "°C", "°C"),
      nc.dimnames = list(x = "lon", y = "lat", time = "time"),
      change.settings = paste0(nc_path, "change_factors.xlsx"),
      save.scenario.matrix = FALSE
  )

}


