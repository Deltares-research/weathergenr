
### Stress test workflow (LOOPING IN SNAKE, NO LONGER IN R)

################################################################################
################################################################################
################################################################################

### INPUTS - EVERYTHING IN THIS SECTION CAN BE PASSED FROM SNAKE
# (DO NOT USE IT IN THE ACTUAL WORKFLOW)

### GENERAL INPUTS (Does not change per run in the loop)

library(weathergenr)

# path to the base nc file [string]
output_path <- "C:/testrun2/future/"

# What prefix we want to attach to the final scenario name when writing back to nc? [string]
nc_file_prefix <- "gabon_climate_runs"

# temp_change_type [string]
temp_change_type = "transient"

# precip_change_type [string]
precip_change_type = "transient"

### RUN SPECIFIC INPUTS (changes in the loop)

# stochastic_nc = Name of the gridded historical realization nc file [string]
stochastic_nc <- "hist_rlz_1.nc"

# vector of monthly precip mean change factors [numeric vector with 12 values]
current_precip_mean_change <- c(0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,   0.7,   0.7,   0.7)

# strtest_matrix_precip_variance = vector of monthly precip variance change factors [numeric vector with 12 values]
current_precip_variance_change <- c(1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,   1.0,   1.0,   1.0)

# strtest_matrix_temp_mean = vector of monthly temperature mean changes, ad DegC [numeric vector with 12 values]
current_temp_mean_change <- c(2.0,  2.0,  2.0,  2.0,  2.0,  2.0,  2.0,  2.0,  2.0,   2.0,  2.0,  2.0)

# What suffix we want to attach to the final scenario name when writing back to nc? [string]
# (index to keep track of both natural variability realization and current climate change run)
nc_file_suffix <- "1_1"

################################################################################
################################################################################
################################################################################


# THIS IS THE MAIN WORKFLOW TO BE CALLED FROM R #

rlz_input_name <- paste0(output_path, "/", stochastic_nc)
rlz_input <- readNetcdf(rlz_input_name)

# Apply climate changes to baseline weather data stored in the nc file
rlz_future <- imposeClimateChanges(
   climate.data = rlz_input$data,
   climate.grid = rlz_input$grid,
   sim.dates = rlz_input$date,
   change.factor.precip.mean = current_precip_mean_change,
   change.factor.precip.variance = current_precip_variance_change,
   change.factor.temp.mean = current_temp_mean_change,
   change.type.temp = temp_change_type,
   change.type.precip = precip_change_type)

 # Save to netcdf file
 writeNetcdf(
   data = rlz_future,
   coord.grid = rlz_input$grid,
   output.path = output_path,
   origin.date =  rlz_input$date,
   calendar.type = "noleap",
   nc.template.file = rlz_input_name,
   nc.compression = 4,
   nc.spatial.ref = "spatial_ref",
   nc.file.prefix = nc_file_prefix,
   nc.file.suffix = nc_file_suffix)


################################################################################


