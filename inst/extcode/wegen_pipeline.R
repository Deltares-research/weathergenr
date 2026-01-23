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
library(weathergenr)
library(dplyr)
library(tidyr)

# =============================================================================
# PROJECT SETUP
# =============================================================================

year_start_month <- 1

# --- Path Configuration ---
#ncfile_dir <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
#out_dir <- file.path("C:/TEMP/ntoum/", year_start_month)
#ncdata <- read_netcdf(nc_path  = paste0(ncfile_dir))

### ---> FOR VECHTE DATA
ncfile_dir   <- "C:/Users/taner/WS/Spongeworks/data/meteo/eobs_v31_1950_2024_allvars_clean.nc"
out_dir <- file.path("C:/TEMP/vechte", year_start_month)
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
  # Variables to simulate
  vars  = c("precip", "temp", "temp_min", "temp_max"),

  # Calendar settings
  year_start_month = year_start_month,    # Water-year start month (1=Jan ... 12=Dec)
  start_year       = 2020,                # First simulation year
  n_years          = NULL,                # Simulated length in years (NULL = all complete years)

  # WARM parameters (annual)
  warm_var         = "precip",            # Variable used for annual wavelet/WARM analysis
  warm_signif      = 0.80,                # Wavelet significance level (for peak detection / filtering)
  warm_pool_size   = 20000L,              # Candidate WARM realizations before filtering

  # Resampling parameters
  n_realizations   = 3L,                  # Number of synthetic realizations to generate
  annual_knn_n     = 100L,                # K for annual KNN analogue selection
  wet_q            = 0.20,                # Wet-day threshold quantile
  extreme_q        = 0.80,                # Extreme-day threshold quantile

  # Spell adjustment factors
  dry_spell_factor = rep(1, 12),          # Multiplier for dry-spell lengths by month
  wet_spell_factor = rep(1, 12),          # Multiplier for wet-spell lengths by month

  # Execution settings
  parallel         = FALSE,               # Enable parallel computation where supported
  n_cores          = NULL,                # Number of cores (NULL = internal default)
  seed             = 100L,                # Base RNG seed (reproducibility)
  verbose          = TRUE,                # Verbose logging
  save_plots       = TRUE                 # Save diagnostic plots to disk (if implemented)
)

config$warm_filter_bounds <- list(
  # --- distributional tolerances (relative diff) ---
  mean                   = 0.03,         # Max |(sim-obs)/obs| for mean
  sd                     = 0.03,         # Max |(sim-obs)/obs| for sd

  # --- tail behaviour (quantile-defined tails + log-distance tol) ---
  tail_low_p             = 0.20,         # Lower-tail quantile (obs threshold)
  tail_high_p            = 0.80,         # Upper-tail quantile (obs threshold)
  tail_tol_log           = log(1.03),    # Max |log(M_sim) - log(M_obs)| for tail mass
  tail_eps               = 1e-5,         # Epsilon for tail-mass log stability

  # --- spectral matching (overall shape) ---
  spectral_cor_min       = 0.60,         # Min correlation of log(GWS_sim) vs log(GWS_obs)

  # --- peak matching (significant observed peaks only) ---
  peak_match_frac_min    = 1.0,          # Minimum fraction of significant obs peaks matched
  n_sig_peaks_max        = 2L,           # Max number of significant obs peaks to enforce
  peak_period_tol        = 0.50,         # Peak period tolerance in log2(period) (octaves)
  peak_mag_tol_log       = log(1.40),    # Max |log(sim/obs)| at matched peaks (magnitude tol)

  # --- plotting diagnostics ---
  plot_wavelet_q         = c(0.05, 0.995), # Quantiles to summarize sim spectra in plots

  # --- relaxation controls ---
  relax_mult             = 1.20,         # Multiplicative relaxation factor for mean/sd/tail tol
  relax_mean_max         = 0.20,         # Upper bound for mean tolerance during relaxation
  relax_sd_max           = 0.20,         # Upper bound for sd tolerance during relaxation
  relax_tail_tol_log_max = log(2.0),     # Upper bound for tail_tol_log during relaxation
  relax_tail_p_step      = 0.02,         # Step for relaxing tail quantile thresholds
  relax_tail_p_low_max   = 0.30,         # Max tail_low_p allowed during relaxation
  relax_tail_p_high_min  = 0.30,         # Min tail_high_p allowed during relaxation

  relax_spectral_cor_step    = 0.05,     # Step decrease for spectral_cor_min
  relax_spectral_cor_min     = 0.30,     # Minimum spectral_cor_min allowed
  relax_peak_match_frac_step = 0.05,     # Step decrease for peak_match_frac_min
  relax_peak_match_frac_min  = 0.00,     # Minimum peak_match_frac_min allowed
  relax_max_iter             = 100L      # Maximum relaxation iterations
)


# =============================================================================
# Step 2: Generate Synthetic Weather
# =============================================================================


stochastic_weather <- generate_weather(
  obs_data           = ncdata$data,
  obs_grid           = ncdata$grid,
  obs_dates          = ncdata$date,
  vars               = config$vars,
  n_years            = config$n_years,
  start_year         = config$start_year,
  year_start_month   = config$year_start_month,
  n_realizations     = config$n_realizations,
  warm_var           = config$warm_var,
  warm_signif        = config$warm_signif,
  warm_pool_size     = config$warm_pool_size,
  warm_filter_bounds = config$warm_filter_bounds,
  annual_knn_n       = config$annual_knn_n,
  wet_q              = config$wet_q,
  extreme_q          = config$extreme_q,
  dry_spell_factor   = config$dry_spell_factor,
  wet_spell_factor   = config$wet_spell_factor,
  out_dir            = out_dir,
  parallel           = config$parallel,
  n_cores            = config$n_cores,
  seed               = config$seed,
  verbose            = config$verbose,
  save_plots         = config$save_plots
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
  daily_sim       = eval_data$sim_data,
  daily_obs       = eval_data$obs_data,
  vars            = config$vars,
  variable_labels = NULL,
  n_realizations  = config$n_realizations,
  eval_max_grids  = 25,
  wet_q           = config$wet_q,
  extreme_q       = config$extreme_q,
  output_dir      = out_dir,
  save_plots      = TRUE,
  seed            = NULL
)



# =============================================================================
# CONVENIENCE WRAPPER
# =============================================================================

out <- run_weather_generator(
  obs_data   = ncdata$data,
  obs_grid   = ncdata$grid,
  obs_dates  = ncdata$date,
  out_dir    = out_dir,
  config     = config,
  eval_max_grids = 25L,
  log_messages = TRUE
)



