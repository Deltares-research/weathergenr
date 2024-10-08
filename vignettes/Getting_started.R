## ----initialize, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE, cache = TRUE)

## ----setup, eval=TRUE, echo = FALSE-------------------------------------------
#devtools::install_github("Deltares/weathergenr", upgrade = "never")
library(weathergenr)

# Other R packages that needs to be installed:
library(dplyr)
library(ggplot2)

## ----basin-image, echo=FALSE, fig.cap="Ntoum basin near Libreville (Gabon)", out.width = '100%'----
knitr::include_graphics('ntoum_basin.png')

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

## ----set_params---------------------------------------------------------------

# Set path to store weather generator results
output_path <- tempdir()

# Set the variables to include in the weather generator
variables <- c("precip", "temp", "temp_min", "temp_max")
variable_labels <- c("Precip.", "Temp. (avg)", "Temp. (min)", "Temp. (max)")
  
# Set the number of stochastic (random) realizations to be generated from the historical weather record
realization_num <- 3

## ----run_wegen, results='hide', eval=TRUE, cache=TRUE-------------------------
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
  warm.sample.num = 5000,
  warm.subset.criteria = NULL,
  knn.sample.num = 100,
  mc.wet.quantile= 0.2,
  mc.extreme.quantile = 0.8,
  output.path = output_path,
  seed = 123)


## ----stochastic3, eval = TRUE, cache=TRUE-------------------------------------
# Resampled dates
stochastic_weather$resampled

# Date vector
head(stochastic_weather$dates)
tail(stochastic_weather$dates)

## ----eval_outputs, eval = TRUE, cache=TRUE, out.width = '100%'----------------

# Find date order from each stochastic simulation
day_order <- sapply(1:realization_num, function(n) match(stochastic_weather$resampled[[n]], ncdata$date))

rlz_sample <- list()
for (n in 1:realization_num) {
  rlz_sample[[n]] <- lapply(ncdata$data[ncdata$grid$id], function(x) x[day_order[,n],] %>%
                              select(precip,  temp, temp_min, temp_max) %>%
                              mutate(date = stochastic_weather$dates, .before = 1))
}

obs_sample <- lapply(ncdata$data[ncdata$grid$id], function(x) x %>%
                       select(precip, temp, temp_min, temp_max) %>%
                       dplyr::mutate(date = ncdata$date, .before = 1))

out <- evaluateWegen(daily.sim = rlz_sample,
              daily.obs = obs_sample,
              output.path = paste0(output_path,"/plots"),
              variables = variables,
              variable.labels = variable_labels,
              variable.units = NULL,
              realization.num = realization_num,
              wet.quantile = 0.2,
              extreme.quantile =0.8)

# Daily means
out$daily_means

# Daily Standard deviations
out$daily_sd

