## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE)

## ----setup, eval=TRUE, echo = FALSE-------------------------------------------
#devtools::install_github("Deltares/weathergenr")
library(weathergenr)

## ----ncfile-------------------------------------------------------------------
ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
ncdata <- readNetcdf(ncfile)

## ----ncdata1------------------------------------------------------------------
names(ncdata)

## ----ncdata2------------------------------------------------------------------
# Display climate data for the first gridcell
ncdata$data[[1]]

## ----ncdata3------------------------------------------------------------------
# Display grid information
ncdata$grid

## ----ncdata4------------------------------------------------------------------
# Display start and ending values date vector
head(ncdata$date)
tail(ncdata$date)

## ----stochastic1--------------------------------------------------------------
# Set path to store weather generator results
output_path <- "C:/testnew/"
variables <- c("precip", "temp", "temp_min", "temp_max")
realization_num <- 1

## ----stochastic2, results='hide', eval = TRUE, cache=TRUE---------------------
start_time <- Sys.time()
stochastic_weather <- generateWeatherSeries(
     weather.data = ncdata$data,
     weather.grid = ncdata$grid,
     weather.date = ncdata$date,
     variable.names = variables,
     variable.labels = variables,
     variable.units = NULL,
     sim.year.num = 20,
     sim.year.start = 2020,
     month.start = 1,
     realization.num = realization_num,
     warm.variable = "precip",
     warm.signif.level = 0.90,
     warm.sample.num = 10000,
     warm.subset.criteria = NULL,
     knn.sample.num = 120,
     mc.wet.quantile= 0.2,
     mc.extreme.quantile = 0.8,
     evaluate.model = TRUE,
     evaluate.grid.num = 20,
     output.path = output_path,
     seed = 1238)
Sys.time() - start_time
## ----stochastic3, eval = TRUE, cache=TRUE-------------------------------------
# Resampled dates
stochastic_weather$resampled

# Date vector
head(stochastic_weather$dates)
tail(stochastic_weather$dates)

## ----climchange, eval = FALSE-------------------------------------------------

 # Temperature changes Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
 delta_temp_mean <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)

 # Precipitation changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
 delta_precip_mean     <- c(0.7, 0.7, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
 delta_precip_variance <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

 # Select first realization
 day_order <- match(stochastic_weather$resampled[[1]], ncdata$date)

 # Obtain stochastic series by re-ordering historical data
 stochastic_rlz <- lapply(ncdata$data, function(x) x[day_order,])

 # Apply climate changes to climate data
 stochastic2 <- imposeClimateChanges(
    climate.data = stochastic_rlz,
    climate.grid = ncdata$grid,
    sim.dates = stochastic_weather$dates,
    change.factor.precip.mean = delta_precip_mean,
    change.factor.precip.variance = delta_precip_variance,
    change.factor.temp.mean = delta_temp_mean,
    transient.temp.change = TRUE,
    transient.precip.change = TRUE,
    calculate.pet = TRUE,
    compute.parallel = FALSE,
    num.cores = NULL,
    fit.method = "mme")

## ----climchange2, eval = FALSE------------------------------------------------
 # Save to netcdf file
 writeNetcdf(
     data = stochastic2,
     coord.grid = ncdata$grid,
     output.path = output_path,
     origin.date =  stochastic_weather$dates[1],
     calendar.type = "noleap",
     nc.template.file = ncfile,
     nc.compression = 4,
     nc.spatial.ref = "spatial_ref",
     nc.file.prefix = "clim",
     nc.file.suffix = NULL)

## ----deltafactors1, eval = FALSE----------------------------------------------
#
 # Temp mean changes Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
 delta_temp_mean_min <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
 delta_temp_mean_max <- c(3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0)

 # Precip mean changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
 delta_precip_mean_min <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
 delta_precip_mean_max <- c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)

 # Precip variance changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
 delta_precip_variance_min <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
 delta_precip_variance_max <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

 # Number of incremental step changes for precip and temp variables
 precip_step_num <- 3
 temp_step_num <- 2

## ----deltafactors2, eval = FALSE----------------------------------------------
 precip_mean_steps <- sapply(1:12, function(m)
          seq(delta_precip_mean_min[m], delta_precip_mean_max[m],
              length.out = precip_step_num))

 precip_variance_steps <- sapply(1:12, function(m)
          seq(delta_precip_variance_min[m], delta_precip_variance_max[m],
              length.out = precip_step_num))

 temp_mean_steps <- sapply(1:12, function(m)
          seq(delta_temp_mean_min[m], delta_temp_mean_max[m],
              length.out = temp_step_num))

  df1 <- as.data.frame(precip_mean_steps) %>% mutate(level = 1:n(),
    variable = "precip_mean", .before = 1)
  df2 <- as.data.frame(precip_variance_steps) %>% mutate(level = 1:n(),
    variable = "precip_variance", .before = 1)
  df3 <- as.data.frame(temp_mean_steps) %>% mutate(level = 1:n(),
    variable = "temp_mean", .before = 1)
  df <- bind_rows(df1, df2, df3) %>% gather(month, value, V1:V12) %>%
    mutate(month = factor(month, levels = paste0("V",1:12), labels = 1:12))

  p <- ggplot2::ggplot(df, aes(x = month, y = value, group = level, color = level)) +
    facet_wrap(. ~ variable, scales = "free_y", ncol = 2) +
    geom_line() +
    labs(x="month", y = "delta factor") +
    scale_color_distiller(palette = "Set1") +
    guides(color = "none")

  p


## ----deltafactors3, eval = FALSE----------------------------------------------
#  # Stress test matrix
strtest_matrix <- tidyr::expand_grid(stoc_ind = 1:realization_num,
  precip_ind = 1:precip_step_num, temp_ind = 1:temp_step_num)

# Total number of scenarios
smax <- nrow(strtest_matrix)

# Stress test delta factors for each variable/climate statistic
strtest_matrix_precip_mean <- precip_mean_steps[strtest_matrix$precip_ind, ]
strtest_matrix_precip_variance <- precip_variance_steps[strtest_matrix$precip_ind, ]
strtest_matrix_temp_mean <- temp_mean_steps[strtest_matrix$temp_ind, ]

## ----deltafactors4, eval = FALSE----------------------------------------------
  write.csv(strtest_matrix,
    paste0(output_path, "strtest_matrix.csv"), row.names = FALSE)
  write.csv(strtest_matrix_precip_mean,
    paste0(output_path, "strtest_matrix_precip_mean.csv"), row.names = FALSE)
  write.csv(strtest_matrix_precip_variance,
    paste0(output_path, "strtest_matrix_precip_variance.csv"), row.names = FALSE)
  write.csv(strtest_matrix_temp_mean,
    paste0(output_path, "strtest_matrix_temp_mean.csv"), row.names = FALSE)

## ----deltafactors5, eval = FALSE----------------------------------------------
  # Read-in resampled dates & date series (from csv files included with the package)
  resampled_dates <- read.csv(system.file("extdata", "resampled_dates.csv", package = "weathergenr"),
    colClasses = "Date")
  sim_dates <- read.csv(system.file("extdata", "sim_dates.csv", package = "weathergenr"),
    colClasses = "Date")[[1]]

  # Use results from generateWeatherSeries function output
  # resampled_dates <- stochastic_weather$resampled
  # sim_dates <- stochastic_weather$dates

 # progress bar (optional)
 pb = txtProgressBar(min = 1, max = smax, initial = 0, style = 3)
 for (s in 1:smax) {

    setTxtProgressBar(pb,s)

    # Find the current scenario indices for the stochastic realization and delta factors
    stoc_ind <- strtest_matrix$stoc_ind[s]

    # Obtain stochastic series by re-ordering historical data
    day_order <- match(resampled_dates[[stoc_ind]], ncdata$date)
    rlz_historical <- lapply(ncdata$data, function(x) x[day_order,])

    # Apply climate changes to climate data
    rlz_future <- imposeClimateChanges(
        climate.data = rlz_historical,
        climate.grid = ncdata$grid,
        sim.dates = sim_dates,
        change.factor.precip.mean = strtest_matrix_precip_mean[s,],
        change.factor.precip.variance = strtest_matrix_precip_variance[s,],
        change.factor.temp.mean = strtest_matrix_temp_mean[s,],
        transient.temp.change = TRUE,
        transient.precip.change = TRUE,
        calculate.pet = TRUE,
        compute.parallel = FALSE,
        num.cores = NULL,
        fit.method = "mme")

      # Save to netcdf file
      # writeNetcdf(
      #   data = rlz_future,
      #   coord.grid = ncdata$grid,
      #   output.path = output_path,
      #   origin.date =  stochastic_weather$dates[1],
      #   calendar.type = "noleap",
      #   nc.template.file = ncfile,
      #   nc.compression = 4,
      #   nc.spatial.ref = "spatial_ref",
      #   nc.file.prefix = "climx",
      #   nc.file.suffix = s)
  }
 close(pb)

