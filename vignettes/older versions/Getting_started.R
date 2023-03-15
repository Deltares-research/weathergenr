## ---- include = FALSE---------------------------------------------------------

library(formatR)

knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE,
  tidy.opts = list(width.cutoff = 60), tidy = TRUE)

## ----setup, eval=TRUE, echo = FALSE-------------------------------------------
#devtools::install_github("Deltares/weathergenr", upgrade = "never")
library(weathergenr)

#Other packages required for the vignette
library(tidyr)
library(dplyr)
library(ggplot2)

## ----basin-image, echo=FALSE, fig.cap="Ntoum basin near Libreville (Gabon)", out.width = '100%'----
knitr::include_graphics('/images/ntoum_basin.png')

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

## ----ncdata5, fig.height = 5, fig.width = 7-----------------------------------
weather_series_basinavg <- (Reduce(`+`, ncdata$data) / nrow(ncdata$grid)) %>%
  as_tibble() %>%
  mutate(date = ncdata$date) %>%
  select(date, precip, temp, temp_min, temp_max) %>%
  gather(key = variable, value = value, -date)

p1 <- ggplot(weather_series_basinavg, aes(date, value, color = variable)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y",ncol = 1) +
  labs(x="Month",y="value") + guides(color = "none")
print(p1)

p2 <- ggplot(weather_series_basinavg,
      aes(x = as.factor(as.numeric(format(date,"%m"))), value, color = variable )) +
  geom_boxplot(aes(group=as.numeric(format(date,"%m")))) +
  facet_wrap(~variable, scales = "free_y",ncol = 2) +
  labs(x="Month",y="value") + guides(color = "none")
print(p2)

## ----stochastic1--------------------------------------------------------------
# Set path to store weather generator results
output_path <- "C:/testrun/intro/"
variables <- c("precip", "temp", "temp_min", "temp_max")
realization_num <- 1

## ----stochastic2, results='hide', eval = FALSE, cache=TRUE--------------------
#  stochastic_weather <- generateWeatherSeries(
#       weather.data = ncdata$data,
#       weather.grid = ncdata$grid,
#       weather.date = ncdata$date,
#       variable.names = variables,
#       variable.labels = variables,
#       variable.units = NULL,
#       sim.year.num = 20,
#       sim.year.start = 2020,
#       month.start = 1,
#       realization.num = realization_num,
#       warm.variable = "precip",
#       warm.signif.level = 0.90,
#       warm.sample.num = 10000,
#       warm.subset.criteria = NULL,
#       knn.sample.num = 120,
#       mc.wet.quantile= 0.2,
#       mc.extreme.quantile = 0.8,
#       evaluate.model = FALSE,
#       evaluate.grid.num = 20,
#       output.path = output_path,
#       seed = 123)

## ----stochastic3, eval = FALSE, cache=TRUE------------------------------------
#  # Resampled dates
#  stochastic_weather$resampled
#
#  # Date vector
#  head(stochastic_weather$dates)
#  tail(stochastic_weather$dates)

## ----climchange, eval = FALSE-------------------------------------------------
#  # Temperature changes Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  delta_temp_mean <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)
#
#  # Precipitation changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  delta_precip_mean     <- c(0.7, 0.7, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#  delta_precip_variance <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
#
#  # Obtain stochastic series by re-ordering historical data
#  day_order <- match(stochastic_weather$resampled[[1]], ncdata$date)
#  stochastic_rlz <- lapply(ncdata$data, function(x) x[day_order,])
#
#  # Apply climate changes to climate data
#  stochastic2 <- imposeClimateChanges(
#     climate.data = stochastic_rlz,
#     climate.grid = ncdata$grid,
#     sim.dates = stochastic_weather$dates,
#     change.factor.precip.mean = delta_precip_mean,
#     change.factor.precip.variance = delta_precip_variance,
#     change.factor.temp.mean = delta_temp_mean,
#     transient.temp.change = TRUE,
#     transient.precip.change = TRUE,
#     calculate.pet = TRUE,
#     compute.parallel = TRUE,
#     num.cores = NULL,
#     fit.method = "mme")

## ----climchange2, eval = FALSE------------------------------------------------
#  # Save to netcdf file
#  writeNetcdf(
#      data = stochastic2,
#      coord.grid = ncdata$grid,
#      output.path = output_path,
#      origin.date =  stochastic_weather$dates[1],
#      calendar.type = "noleap",
#      nc.template.file = ncfile,
#      nc.compression = 4,
#      nc.spatial.ref = "spatial_ref",
#      nc.file.prefix = "clim",
#      nc.file.suffix = NULL)

