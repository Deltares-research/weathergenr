
#load("./ntoum_m1.Rdata")

# Libraries
library(devtools)
library(weathergenr)
library(dplyr)
library(tidyr)
library(microbenchmark)


# Read-in forcing data & geometry
month.start <- 1
output_path <- paste0("C:/TEMP/", month.start, "/")


#case_path <- "C:/Users/taner/WS/SpongeWorks/"
#ncfile <- paste0(case_path, "data/meteo/extract_historical.nc")
#ncdata <- read_netcdf(ncfile)

#case_path <- "C:/Users/taner/WS/SpongeWorks/"
#nc_eops <- paste0(case_path, "data/meteo/eobs_v31_1950_2024_allvars_clean.nc")
#nc_file <- nc_eops
#ncdata <- read_netcdf(nc_file, variables = c("precip", "temp", "tn", "tx"),
#                      var_rename = c(tn = "temp_min", tx = "temp_max"))

# Test using default ntoum data (Liberia)
ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
ncdata <- read_netcdf(ncfile)


#### Define all variables in advance for testing
mc.wet.quantile <- 0.1
weather.data <- ncdata$data
weather.grid <- ncdata$grid
weather.date <- ncdata$date
variables <- c("precip", "temp", "temp_min", "temp_max")
variable.labels <- variables
sim.year.num <- 20
sim.year.start <- 2020
realization.num <- 3
warm.variable <- "precip"
warm.signif.level <- 0.80
warm.sample.num <- 10000
warm.subset.criteria <- list(mean = 0.05, sd = 0.05, min = 0.05, max = 0.05, sig.thr = 0.8, nsig.thr = 1.5)
knn.sample.num <- 100
mc.extreme.quantile <- 0.8
dry.spell.change <- rep(1, 12)
wet.spell.change <- rep(1, 12)
output.path <- output_path
seed <- 1000
compute.parallel<- FALSE
num.cores <- NULL


################################################################################

stochastic_weather <- generateWeatherSeries(
  weather.data = weather.data,
  weather.grid = weather.grid,
  weather.date = weather.date,
  variables = variables,
  sim.year.num = sim.year.num,
  sim.year.start = sim.year.start,
  month.start = month.start,
  realization.num = realization.num,
  warm.variable = warm.variable,
  warm.signif.level = warm.signif.level,
  warm.sample.num = warm.sample.num,
  warm.subset.criteria = warm.subset.criteria,
  knn.sample.num = knn.sample.num,
  mc.wet.quantile= mc.wet.quantile,
  mc.extreme.quantile = mc.extreme.quantile,
  dry.spell.change = dry.spell.change,
  wet.spell.change = wet.spell.change,
  output.path = output.path,
  compute.parallel = FALSE,
  num.cores = num.cores,
  seed = seed)

day_order <- sapply(1:realization.num,
     function(n) match(stochastic_weather$resampled[[n]], ncdata$date))


sim_date_complete_years <- tibble(date = stochastic_weather$dates) %>%
  mutate(year = as.integer(format(date, "%Y"))) %>%
  group_by(year) %>%
  filter(n() >= 365) %>%
  ungroup()

sim_date_complete_years_idx <- which(stochastic_weather$dates %in% sim_date_complete_years$date)

rlz_sample <- list()
for (n in 1:realization.num) {
  rlz_sample[[n]] <- lapply(ncdata$data[ncdata$grid$id], function(x) x[day_order[,n],] %>%
                              select(precip,  temp, temp_min, temp_max) %>%
                              mutate(date = stochastic_weather$dates, .before = 1) %>%
                              slice(sim_date_complete_years_idx))
}

obs_sample <- lapply(ncdata$data[ncdata$grid$id], function(x) x %>%
                       select(precip, temp, temp_min, temp_max) %>%
                       dplyr::mutate(date = ncdata$date, .before = 1))


out <- evaluate_weather_generator(daily.sim = rlz_sample,
                     daily.obs = obs_sample,
                     output.path = output_path,
                     save.plots = TRUE,
                     variables = variables,
                     variable.labels = variables,
                     variable.units = NULL,
                     realization.num = realization.num,
                     wet.quantile = mc.wet.quantile,
                     extreme.quantile = mc.extreme.quantile)
