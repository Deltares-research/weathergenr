# =============================================================================
# Stochastic Weather Generation Workflow
# =============================================================================
#
# This script demonstrates the complete workflow for:
#   1. Reading climate data from NetCDF files
#   2. Generating synthetic weather using the WARM methodology
#   3. Evaluating generator performance against observations
#
# =============================================================================

# Load required libraries
library(devtools)
library(weathergenr)
library(dplyr)
library(tidyr)

# =============================================================================
# PROJECT SETUP
# =============================================================================

year_start_month <- 10

# --- Path Configuration ---
#ncfile_dir <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
#output_dir <- file.path("C:/TEMP/ntoum/", year_start_month)
#ncdata <- read_netcdf(nc_path  = paste0(ncfile_dir))

### ---> FOR VECHTE DATA
ncfile_dir   <- "C:/Users/taner/WS/Spongeworks/data/meteo/eobs_v31_1950_2024_allvars_clean.nc"
output_dir <- file.path("C:/TEMP/vechte/", year_start_month)
ncdata <- read_netcdf(nc_path  = paste0(ncfile_dir),
    var = c("precip", "temp", "tn", "tx"),
    var_name = c(tn = "temp_min", tx = "temp_max"))


### ---> FOR RHINE BASIN
#ncfile_dir   <- "C:/Users/taner/WS/Spongeworks/data/meteo/eobs_v31_1950_2024_allvars_clean.nc"
#output_dir <- "C:/TEMP/vechte/"

# =============================================================================
# KEY PARAMETERS
# =============================================================================

# --- Key Simulation Parameters ---
config <- list(
  # Calendar settings
  year_start_month = year_start_month,   # Water year starts in October (1 = calendar year)
  start_year       = 2020,       # First simulation year
  n_years          = NULL,       # NULL = use all available complete years

  # Variables to simulate
  vars = c("precip", "temp", "temp_min", "temp_max"),

  # WARM parameters
  warm_var       = "precip",     # Variable for annual wavelet analysis
  warm_signif    = 0.80,         # Wavelet significance threshold
  warm_pool_size = 20000,        # Candidate realizations before filtering

  # Resampling parameters
  n_realizations = 3,            # Number of synthetic series to generate
  annual_knn_n   = 100,          # K for annual KNN matching
  wet_q          = 0.20,         # Wet day threshold quantile
  extreme_q      = 0.80,         # Extreme day threshold quantile

  # Spell adjustment factors (length 12, one per month)
  dry_spell_factor = rep(1, 12),
  wet_spell_factor = rep(1, 12),

  # Execution settings
  parallel = FALSE,
  n_cores  = NULL,
  seed     = 1030,
  verbose  = TRUE
)

# =============================================================================
# Step 2: Generate Synthetic Weather
# =============================================================================

stochastic_weather <- generate_weather(
  obs_data         = ncdata$data,
  obs_grid         = ncdata$grid,
  obs_dates        = ncdata$date,
  vars             = config$vars,
  n_years          = config$n_years,
  start_year       = config$start_year,
  year_start_month = config$year_start_month,
  n_realizations   = config$n_realizations,
  warm_var         = config$warm_var,
  warm_signif      = config$warm_signif,
  warm_pool_size   = config$warm_pool_size,
  annual_knn_n     = config$annual_knn_n,
  wet_q            = config$wet_q,
  extreme_q        = config$extreme_q,
  dry_spell_factor = config$dry_spell_factor,
  wet_spell_factor = config$wet_spell_factor,
  out_dir          = output_dir,
  parallel         = config$parallel,
  n_cores          = config$n_cores,
  seed             = config$seed,
  verbose          = config$verbose
)

# =============================================================================
# Step 3: Prepare Data for Evaluation
# =============================================================================

eval_data <- prepare_evaluation_data(
  gen_output = stochastic_weather,
  obs_data   = ncdata$data,
  obs_dates  = ncdata$date,
  grid_ids   = ncdata$grid$id,
  variables  = config$vars,
  verbose = TRUE
)

# =============================================================================
# Step 4: Evaluate Generator Performance
# =============================================================================

evaluation <- evaluate_weather_generator(
  daily_sim        = eval_data$sim_data,
  daily_obs        = eval_data$obs_data,
  variables        = config$vars,
  variable_labels  = NULL,
  n_realizations   = config$n_realizations,
  wet_quantile     = config$wet_q,
  extreme_quantile = config$extreme_q,
  output_dir      = output_dir,
  save_plots       = TRUE,
  max_grids        = 25,
  seed             = NULL
)
