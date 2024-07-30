

library(dplyr)
library(ggplot2)
library(tidyr)
library(weathergenr)

output_path <- "C:\\Users\\taner\\Workspace\\test-weathergenr\\data\\vakhsh\\"
ncfile <- paste0(output_path, "extract_historical.nc")
ncdata <- readNetcdf(ncfile)

## ----stochastic1--------------------------------------------------------------

# Set path to store weather generator results
variables <- c("precip", "temp", "temp_min", "temp_max")
realization_num <- 15

## ----stochastic2, results='hide', eval = TRUE, cache=TRUE---------------------

warm_criteria = list(
   mean = c(0.90,1.10),
   sd = c(0.85,1.15),
   min = c(0.90,1.10),
   max = c(0.90,1.10),
   power = c(0.50,10),
   nonsignif.threshold = 1.20)

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
  warm.signif.level = 0.80,
  warm.sample.num = 30000,
  warm.subset.criteria = warm_criteria,
  knn.sample.num = 100,
  mc.wet.quantile= 0.2,
  mc.extreme.quantile = 0.8,
  evaluate.model = TRUE,
  evaluate.grid.num = 30,
  output.path = output_path,
  compute.parallel = TRUE,
  seed = 555)



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
