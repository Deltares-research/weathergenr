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
#out_dir <- file.path("C:/TEMP/ntoum/", year_start_month)
#ncdata <- read_netcdf(nc_path  = paste0(ncfile_dir))

### ---> FOR VECHTE DATA
ncfile_dir   <- "C:/Users/taner/WS/Spongeworks/data/meteo/eobs_v31_1950_2024_allvars_clean.nc"
out_dir <- file.path("C:/TEMP/vechte/", year_start_month)
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
  start_year       = 2020,               # First simulation year
  n_years          = NULL,               # NULL = use all available complete years

  # Variables to simulate
  vars = c("precip", "temp", "temp_min", "temp_max"),

  # WARM parameters
  warm_var       = "precip",     # Variable for annual wavelet analysis
  warm_signif    = 0.80,         # Wavelet significance threshold
  warm_pool_size = 50000,        # Candidate realizations before filtering
  warm_filter_relax_order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),

  # Resampling parameters
  n_realizations = 3,            # Number of synthetic series to generate
  annual_knn_n   = 100,          # K for annual KNN matching
  wet_q          = 0.20,         # Wet day threshold quantile
  extreme_q      = 0.80,         # Extreme day threshold quantile

  # Spell adjustment factors (length 12, one per month)
  dry_spell_factor = rep(1, 12),
  wet_spell_factor = rep(1, 12),

  # Execution settings
  parallel = TRUE,
  n_cores  = NULL,
  seed     = 1030,
  verbose  = TRUE,

  warm_filter_bounds = list(

    # --- distributional tolerances (relative diff) ---
    mean = 0.03,
    sd   = 0.03,

    # --- tail behaviour (quantile-defined tails + log-distance tol) ---
    tail_low_p   = 0.20,
    tail_high_p  = 0.80,
    tail_tol_log = log(1.03),
    tail_eps     = 1e-5,

    # --- spectral matching (overall shape) ---
    spectral_cor_min = 0.60,   # Min correlation of log-GWS

    # --- peak matching (significant observed peaks only) ---
    peak_match_frac_min = 1.0,        # Fraction of significant observed peaks that must match
    n_sig_peaks_max     = 2,          # Max number of significant observed peaks to enforce
    peak_period_tol     = 0.50,       # Period tolerance in log2 scale (octaves)
    peak_mag_tol_log    = log(1.40),  # abs(log(sim/obs)) <= log(1.5) (~ within +/-50%)

    # --- plotting diagnostics ---
    plot_wavelet_q = c(0.05, 0.995),

    # --- relaxation controls ---
    relax_mult     = 1.20,
    relax_mean_max = 0.20,
    relax_sd_max   = 0.20,

    relax_tail_tol_log_max = log(2.0),
    relax_tail_p_step      = 0.02,
    relax_tail_p_low_max   = 0.30,
    relax_tail_p_high_min  = 0.30,

    # Spectral relaxation (simplified)
    relax_spectral_cor_step    = 0.05,
    relax_spectral_cor_min     = 0.30,
    relax_peak_match_frac_step = 0.05,
    relax_peak_match_frac_min  = 0.00,

    relax_max_iter = 100L
  )
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
  warm_filter_bounds = config$warm_filter_bounds,
  annual_knn_n     = config$annual_knn_n,
  wet_q            = config$wet_q,
  extreme_q        = config$extreme_q,
  dry_spell_factor = config$dry_spell_factor,
  wet_spell_factor = config$wet_spell_factor,
  out_dir          = out_dir,
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
  output_dir      = out_dir,
  save_plots       = TRUE,
  max_grids        = 25,
  seed             = NULL
)
