
library(weathergenr)
library(dplyr)
library(tidyr)
library(ggplot2)

ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
ncdata <- read_netcdf(ncfile, var_name = c("precip" = "prcp"))

prcp_obs <- ncdata$data[[1]]$prcp
obs_date <- ncdata$date
obs_month <- as.integer(format(obs_date, "%m"))
obs_year  <- as.integer(format(obs_date, "%Y"))
obs_year <- obs_year - min(obs_year) + 1
n_years <- max(obs_year)

prcp_mean_change <- rep(0.7, 12)
prcp_var_change <- rep(1.0, 12)
prcp_mean_change_m <- matrix(prcp_mean_change, nrow = n_years, ncol = 12, byrow = TRUE)
prcp_var_change_m <- matrix(prcp_var_change, nrow = n_years, ncol = 12, byrow = TRUE)

prcp_qm <- adjust_precipitation_qm(
  prcp = prcp_obs,
  mean_factor = prcp_mean_change_m,
  var_factor = prcp_var_change_m,
  scale_var_with_mean = TRUE,
  exaggerate_extremes = TRUE,
  enforce_target_mean = TRUE,
  month = obs_month,
  year = obs_year,
  min_events = 10,
  validate_output = TRUE,
  diagnostics = TRUE,
  verbose = FALSE
)








prcp_qm$diagnostics
# Apply perturbations (returns the perturbed realization)
rlz_future <- apply_climate_perturbations(
  data = ncdata$data,
  grid = ncdata$grid,
  date = ncdata$date,
  prcp_mean_factor = prcp_mean_factor_s,
  prcp_var_factor  = prcp_var_factor_s,
  temp_delta       = temp_delta_s,
  temp_transient   = TRUE,
  prcp_transient   = TRUE,
  compute_pet      = FALSE,
  verbose          = FALSE
)













#Keep this small for a quick vignette run.
n_realizations <- 2
set.seed(1)

sim_weather <- generate_weather(
  obs_data = ncdata$data,
  obs_grid = ncdata$grid,
  obs_dates = ncdata$date,
  vars = weather_vars,
  n_years = 20,
  start_year = 2020,
  year_start_month = 1,
  n_realizations = n_realizations,
  warm_var = "precip",
  warm_signif = 0.90,
  warm_pool_size = 5000,
  annual_knn_n = 30,
  wet_q = 0.2,
  extreme_q = 0.8,
  out_dir = out_dir,
  seed = 123,
  parallel = FALSE,
  verbose = FALSE
)


#Monthly ranges (12 values). Replace these with your own stress-test design.
temp_mean_delta_min <- rep(0.0, 12)
temp_mean_delta_max <- rep(3.0, 12)

prcp_mean_factor_min <- rep(0.7, 12)
prcp_mean_factor_max <- rep(1.3, 12)

#Keep variance fixed here; set a range (e.g., 0.8-1.5) if you want to stress extremes.
prcp_var_factor_min <- rep(1.0, 12)
prcp_var_factor_max <- rep(1.0, 12)

#Discretization level (number of scenario levels per axis)
n_prcp_levels <- 3
n_temp_levels <- 2

prcp_mean_factor_levels <- matrix(NA_real_, nrow = n_prcp_levels, ncol = 12)
prcp_var_factor_levels <- matrix(NA_real_, nrow = n_prcp_levels, ncol = 12)
temp_mean_delta_levels <- matrix(NA_real_, nrow = n_temp_levels, ncol = 12)

for (m in 1:12) {
  prcp_mean_factor_levels[, m] <- seq(prcp_mean_factor_min[m], prcp_mean_factor_max[m], length.out = n_prcp_levels)
  prcp_var_factor_levels[, m] <- seq(prcp_var_factor_min[m], prcp_var_factor_max[m], length.out = n_prcp_levels)
  temp_mean_delta_levels[, m] <- seq(temp_mean_delta_min[m], temp_mean_delta_max[m], length.out = n_temp_levels)
}

colnames(prcp_mean_factor_levels) <- paste0("V", 1:12)
colnames(prcp_var_factor_levels) <- paste0("V", 1:12)
colnames(temp_mean_delta_levels) <- paste0("V", 1:12)

df_prcp_mean <- cbind(variable = "prcp_mean", level = seq_len(n_prcp_levels), as.data.frame(prcp_mean_factor_levels))
df_prcp_var <- cbind(variable = "prcp_variance", level = seq_len(n_prcp_levels), as.data.frame(prcp_var_factor_levels))
df_temp_mean <- cbind(variable = "temp_mean", level = seq_len(n_temp_levels), as.data.frame(temp_mean_delta_levels))

scenario_index <- expand.grid(
  stoc_ind = 1:n_realizations,
  precip_ind = 1:n_prcp_levels,
  temp_ind = 1:n_temp_levels
)

scenario_prcp_mean_factor <- prcp_mean_factor_levels[scenario_index$precip_ind, ]
scenario_prcp_var_factor  <- prcp_var_factor_levels[scenario_index$precip_ind, ]
scenario_temp_mean_delta  <- temp_mean_delta_levels[scenario_index$temp_ind, ]



# Example: pick a single scenario row (e.g., the first row)
s <- 1

# Which realization to perturb
r <- scenario_index$stoc_ind[s]

# Monthly perturbations for this scenario row (length-12 vectors)
prcp_mean_factor_s <- as.numeric(scenario_prcp_mean_factor[s, ])
prcp_var_factor_s  <- as.numeric(scenario_prcp_var_factor[s, ])
temp_delta_s       <- as.numeric(scenario_temp_mean_delta[s, ])

# Select one realization from the stochastic output
sim_dates <- sim_weather$dates
sim_dates_resampled <- sim_weather$resampled[[r]]

# Create day order indices for each realization
sim_day_order <- match(sim_dates_resampled, ncdata$date)

# Extract synthetic weather samples
rlz_historical <- lapply(ncdata$data[ncdata$grid$id], function(x) {
    x[sim_day_order, ] %>% select(prcp = precip, temp, temp_min, temp_max) })

prcp_org <- rlz_historical[[1]]$prcp
sim_month <- lubridate::month(sim_dates)
sim_year <- lubridate::year(sim_dates)
sim_year <- sim_year - min(sim_year) + 1

prcp_mean_factor <- matrix(prcp_mean_factor_s, nrow = max(sim_year), ncol = 12)
prcp_var_factor  <- matrix(prcp_var_factor_s, nrow = max(sim_year), ncol = 12)



# Apply perturbations (returns the perturbed realization)
rlz_future <- apply_climate_perturbations(
  data = rlz_historical,
  grid = ncdata$grid,
  date = sim_dates,
  prcp_mean_factor = prcp_mean_factor_s,
  prcp_var_factor  = prcp_var_factor_s,
  temp_delta       = temp_delta_s,
  temp_transient   = TRUE,
  prcp_transient   = TRUE,
  compute_pet      = FALSE,
  verbose          = FALSE
)



