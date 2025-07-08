# ===================================
# Load required packages
# ===================================
library(weathergenr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

# ===================================
# Set global parameters
# ===================================
output_path <- "C:/TEMP"
realization_num <- 3
precip_multipliers <- c(0.75, 1, 1.25)
precip_labels <- paste0("x", precip_multipliers)
variables <- c("precip", "temp", "temp_min", "temp_max")

# ===================================
# Load and read NetCDF data
# ===================================
ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
ncdata <- readNetcdf(ncfile)

# ===================================
# Generate synthetic weather series
# ===================================
res <- generateWeatherSeries(
  weather.data         = ncdata$data,
  weather.grid         = ncdata$grid,
  weather.date         = ncdata$date,
  variable.names       = variables,
  variable.labels      = variables,
  variable.units       = NULL,
  sim.year.num         = 20,
  sim.year.start       = 2020,
  month.start          = 1,
  realization.num      = realization_num,
  warm.variable        = "precip",
  warm.signif.level    = 0.90,
  warm.sample.num      = 10000,
  warm.subset.criteria = NULL,
  knn.sample.num       = 100,
  mc.wet.quantile      = 0.3,
  mc.extreme.quantile  = 0.8,
  output.path          = output_path,
  seed                 = 555,
  compute.parallel     = FALSE,
  num.cores            = NULL,
  save.rdata           = FALSE
)

# Obtain stochastic series by re-ordering historical data
rlz_historical <- list()
for (x in seq_len(realization_num)) {
  day_order <- match(res$resampled[[x]], ncdata$date)
  rlz_historical[[x]] <- lapply(ncdata$data, function(y) y[day_order,])
}

# Create 3 scenarios by multiplying precipitation
precip_multipliers <- c(0.75,1,1.25)
precip_multipliers_labs <- paste0("x",precip_multipliers)


# Loop to apply precipitation changes over each realization
rlz_future <- list()
precip_future <- list()

for (x in seq_len(realization_num)) {

  for (y in seq_len(ncol(mat_change_precip_means))) {

    # Apply climate changes to climate data
    rlz_future[[y]] <- imposeClimateChanges(
      climate.data = rlz_historical[[x]],
      climate.grid = ncdata$grid,
      sim.dates = res$dates,
      change.factor.precip.mean = rep(precip_multipliers[y],12),
      change.factor.precip.variance = rep(1, 12),
      change.factor.temp.mean = rep(0, 12),
      transient.temp.change = TRUE,
      transient.precip.change = FALSE,
      calculate.pet = TRUE,
      compute.parallel = FALSE,
      num.cores = NULL,
      fit.method = "mme")

    # Save to netcdf file
    writeNetcdf(
      data = rlz_future[[y]],
      coord.grid = ncdata$grid,
      output.path = output_path,
      origin.date =  res$dates[1],
      calendar.type = "noleap",
      nc.template.file = ncfile,
      nc.compression = 4,
      nc.spatial.ref = "spatial_ref",
      nc.file.prefix = "climx",
      nc.file.suffix = 1
    )
  }

  # Aggregate to Monthly totals
  df <- lapply(rlz_future, function(z) z[[1]])
  precip_future[[x]] <- do.call(cbind,
                                lapply(df, function(df) df["precip"])) %>%
    setNames(precip_multipliers_labs) %>% as_tibble() %>%
    mutate(date = res$dates) %>%
    pivot_longer(names_to = "scenario", values_to = "value", cols = -date) %>%
    mutate(year = lubridate::year(date), month = lubridate::month(date)) %>%
    group_by(year, month, scenario) %>%
    summarise(value = sum(value, na.rm = TRUE))

}

# ===================================
# Combine all realizations and plot
# ===================================
precip_future_comb <- bind_rows(precip_future, .id = "rlz") %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

my_colors <- c("x0.75" = "red", "x1" = "#000000", "x1.25" = "royalblue")

ggplot(precip_future_comb, aes(x = date, y = value, color = scenario, group = scenario)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~ rlz, nrow = realization_num) +
  scale_color_manual(values = my_colors) +
  labs(
    title = "Monthly Precipitation Under Different Scenarios",
    x = "Date",
    y = "Monthly Precipitation (mm)",
    color = "Scenario"
  ) +
  theme_bw()
