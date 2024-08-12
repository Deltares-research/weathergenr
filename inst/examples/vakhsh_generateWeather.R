
#### Script to generate synthetic weather traces using a gridded weather forcing dataset (netcdf file)
#### Last updated: 9/8/2024

# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(weathergenr)
library(sf)
library(sfheaders)

# Additional settings
sf_use_s2(FALSE)

## ----Specify forcing dataset -------------------------------------------------

data_path <- "C:/Users/taner/Workspace/weathergenr-applications/Vakhsh/20240809/"
output_path <- paste0(data_path, "output/")

ncfile <- paste0(data_path, "extract_historical.nc")
ncdata <- readNetcdf(ncfile)

## ----Set key variables -------------------------------------------------------

# Variables to include
variables <- c("precip", "temp", "temp_min", "temp_max")

# Number of stochastic realizations to generate
realization_num <- 11

# Wavelet AR model filtering criteria
warm_criteria = list(
  mean = c(0.95,1.05),
  sd = c(0.95,1.05),
  min = c(0.90,1.10),
  max = c(0.90,1.10),
  power = c(0.80,10),
  nonsignif.threshold = 1.20)

## ----stochastic2, results='hide', eval = TRUE, cache=TRUE---------------------

stochastic_weather <- generateWeatherSeries(
  weather.data = ncdata$data,
  weather.grid = ncdata$grid,
  weather.date = ncdata$date,
  variable.names = variables,
  variable.labels = variables,
  variable.units = NULL,
  sim.year.num = 52.0,
  sim.year.start = 1980,
  month.start = 1,
  realization.num = realization_num,
  warm.variable = "precip",
  warm.signif.level = 0.95,
  warm.sample.num = 30000,  # suggested range 10,000 - 50,000
  warm.subset.criteria = warm_criteria,
  knn.sample.num = 120, # suggested 100 or 120
  mc.wet.quantile = 0.2, # don't change
  mc.extreme.quantile = 0.8, # don't change
  output.path = output_path,
  compute.parallel = FALSE,
  seed = 2024  # Randomization seed
)

################################################################################
## ----Evaluate results, results='hide', eval = TRUE, cache=TRUE----------------

# Read in historical and synthetic weather data
ncfile <- paste0(data_path, "extract_historical.nc")
sim_dates <- readr::read_csv(paste0(output_path, "sim_dates.csv"))
resampled_dates <- readr::read_csv(paste0(output_path, "resampled_dates.csv"))
ncdata <- weathergenr::readNetcdf(ncfile)
realization.num <- ncol(resampled_dates)
ngrids <- length(ncdata$data)
weather.date <- ncdata$date
weather.data <- ncdata$data
weather.grid <- ncdata$grid

# Read-in geometry
rivers_sf <- read_sf(paste0(data_path, "staticgeoms/rivers.geojson"))
basins_sf <- read_sf(paste0(data_path, "staticgeoms/basins.geojson"))
basins_vakhsh_sf <- read_sf(paste0(data_path, "staticgeoms/basins_vakhsh.geojson"))

evaluate.grid.num <- 50

# Sample evenly from the grid cells using sf package
gridpoints_sf <- sf::st_as_sf(weather.grid, coords=c("x","y"), remove = FALSE) %>%
  sf::st_set_crs(4326)

# Intersecting points
basinpoints <- which(st_intersects(gridpoints_sf, basins_sf, sparse = FALSE))

gridpoints_sf_sample <- gridpoints_sf %>%
  filter(id %in% basinpoints) %>%
  sf::st_sample(size = min(evaluate.grid.num, ngrids), type="regular") %>%
  sf::st_cast("POINT") %>%
  sf::st_coordinates() %>%
  as_tibble() %>%
  rename(c("x"="X","y"="Y")) %>%
  inner_join(gridpoints_sf,., by = c("x","y"))

p <- ggplot() +
  theme_light() +
  geom_sf(data = basins_sf, fill = "lightgreen", color = "black",  alpha = 0.3) +
  geom_sf(data = rivers_sf, color = "black", , alpha = 0.3) +
  geom_sf(data = basins_vakhsh_sf, fill = "blue", alpha = 0.3) +
  geom_sf(data = gridpoints_sf, color = "gray", size = 1, alpha = 0.5) +
  geom_sf(data = gridpoints_sf_sample, color = "red", size = 3, shape = 8)

ggsave(filename = paste0(output_path,"eval_sample_grids.png"), width = 7, height = 5)
sampleGrids <- gridpoints_sf_sample$id

day_order <- sapply(1:realization.num,
      function(n) match(resampled_dates[[n]], weather.date))


rlz_sample <- list()
for (n in 1:realization.num) {
  rlz_sample[[n]] <- lapply(weather.data[sampleGrids], function(x) x[day_order[,n],] %>%
                              select(precip,  temp, temp_min, temp_max) %>%
                              mutate(date = sim_dates[[1]], .before = 1))
}

obs_sample <- lapply(weather.data[sampleGrids], function(x) x %>%
                       select(precip, temp, temp_min, temp_max) %>%
                       dplyr::mutate(date = weather.date, .before = 1))

out <- evaluateWegen(daily.sim = rlz_sample,
              daily.obs = obs_sample,
              output.path = paste0(output_path,"plots1"),
              variables = c("precip", "temp", "temp_min", "temp_max"),
              variable.labels = c("Precip.", "Temp. (avg)", "Temp. (min)", "Temp. (max)"),
              variable.units = NULL,
              realization.num = realization.num,
              wet.quantile = 0.2,
              extreme.quantile =0.8)



################################################################################

## ----clim change, eval = FALSE------------------------------------------------

# Temperature changes Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
delta_temp_mean <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)

# Precipitation changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
delta_precip_mean     <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
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
  transient.temp.change = FALSE,
  transient.precip.change = FALSE,
  calculate.pet = TRUE,
  compute.parallel = FALSE,
  num.cores = 1,
  fit.method = "mme")


# Non-exceedance probability
library(dplyr)
library(tidyr)


df <- tibble(date = stochastic_weather$dates) %>%
  mutate(
    month = as.numeric(format(date,"%m"),
    scn = stochastic2[[1]]$precip,
         base = stochastic_rlz[[1]]$precip))


sim_cc <- stochastic2[[1]] %>%
  mutate(date = stochastic_weather$dates)

sim_base <- stochastic_rlz[[1]] %>%
  mutate(date = stochastic_weather$dates)

sims <- sim_base %>% select(date, no_change = precip) %>%
  mutate(cli_change = sim_cc$precip) %>%

  select(-date)


df <- sims %>%
  gather(key = variable, value = value, -month) %>%
  group_by(variable, month) %>%
  mutate(value = sort(value, decreasing = FALSE, na.last = TRUE),
         x = 100*(1:n())/(n()+1)) %>%
  filter(value > 0)

var_levels <- unique(df$variable)
var_colors <- c("black", "blue", "red", "yellow")[1:length(var_levels)]

#Print ggplot object
p <- ggplot(df, aes(x = x, y = value)) +
  geom_line(aes(color = variable, group = variable)) +
  scale_color_manual(values = var_colors)  +
  facet_wrap(~ month, scales = "free_y") +
  xlab("Percentage exceedance (%)")


###############################################################################
###############################################################################

ggFDC(sims$no_change, sims$cli_change)

#Flow-duration-curves
ggFDC <- function(..., log.scale = FALSE, panel.var = "month") {

  require(dplyr)
  require(ggplot2)
  require(tidyr)

  #Dataframe
  df <- data_frame(...)

  df <- df %>%
    gather(key = variable, value = value) %>%
    group_by(variable, month) %>%
    mutate(value = sort(value, decreasing = TRUE, na.last = TRUE),
           x = 100*(1:n())/(n()+1)) %>%
    filter(value > 0)

  var_levels <- unique(df$variable)
  var_colors <- c("black", "blue", "red", "yellow")[1:length(var_levels)]

  #Print ggplot object
  p <- ggplot(df, aes(x = x, y = value)) +
    geom_line(aes(color = variable, group = variable)) +
    scale_color_manual(values = var_colors)  +
    facet_wrap(~ month, scales = "free_y") +
    xlab("Percentage exceedance (%)") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))


  if(log.scale == TRUE) {

    require(scales)

    p + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))
  }

  return(p)
}


ggFDC(sims$no_change, sims$cli_change, panel.var = "month")

library(dplyr)
library(ggplot2)

weather.data = ncdata$data
weather.grid = ncdata$grid
weather.date = ncdata$date
variable.names = variables
variable.labels = variables
variable.units = NULL
sim.year.num = 30
sim.year.start = 2010
month.start = 1
realization.num = realization_num
warm.variable = "precip"
warm.signif.level = 0.90
warm.sample.num = 10000
warm.subset.criteria = NULL
knn.sample.num = 120
mc.wet.quantile= 0.2
mc.extreme.quantile = 0.8
evaluate.model = FALSE
evaluate.grid.num = 20
output.path = output_path
compute.parallel = FALSE
seed = 1238
