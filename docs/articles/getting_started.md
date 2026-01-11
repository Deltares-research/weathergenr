# Getting Started with weathergenr

## Introduction

The **weathergenr** package provides a comprehensive toolkit for
stochastic weather generation and climate data analysis. It implements
sophisticated statistical methods to create synthetic daily weather
series that preserve the statistical properties of observed climate data
while enabling scenario-based climate impact assessments.

This vignette demonstrates the core workflow:

1.  Reading gridded meteorological data from NetCDF files
2.  Generating synthetic weather realizations using coupled statistical
    methods  
3.  Evaluating the quality of synthetic weather against historical
    observations

### Load Required Packages

Code

``` r

library(weathergenr)
library(dplyr)
library(ggplot2)
```

> **Related functions**
>
> - \[[`read_netcdf()`](https://github.com/Deltares-research/weathergenr/reference/read_netcdf.md)\]
> - \[[`generate_weather_series()`](https://github.com/Deltares-research/weathergenr/reference/generate_weather_series.md)\]
> - \[[`evaluate_weather_generator()`](https://github.com/Deltares-research/weathergenr/reference/evaluate_weather_generator.md)\]
> - \[[`write_netcdf()`](https://github.com/Deltares-research/weathergenr/reference/write_netcdf.md)\]
> - \[[`apply_climate_perturbations()`](https://github.com/Deltares-research/weathergenr/reference/apply_climate_perturbations.md)\]
> - \[[`wavelet_spectral_analysis()`](https://github.com/Deltares-research/weathergenr/reference/wavelet_spectral_analysis.md)\]

## 📂 Reading Gridded Multivariate Weather Data

### About the Example Dataset

For this tutorial, we use the ERA5 dataset (ECMWF Reanalysis v5)
included with the package. The data covers the Ntoum Basin in Gabon, a
tropical watershed in West Central Africa. The dataset includes daily
precipitation, mean temperature, minimum temperature, and maximum
temperature at a 0.25-degree spatial resolution.

![](figures/ntoum_basin.png)

Map of Ntoum Basin, Gabon

### Loading NetCDF Data

The
[`read_netcdf()`](https://github.com/Deltares-research/weathergenr/reference/read_netcdf.md)
function provides an efficient interface for extracting meteorological
data from NetCDF files. It handles the complexities of NetCDF format and
returns data in a tidy, analysis-ready structure optimized for the
weathergenr workflow.

#### Key Features

The function performs several important operations:

- Extracts spatial coordinates and temporal information
- Organizes data into per-grid-cell time series
- Handles calendar systems including leap day management
- Provides variable renaming capabilities
- Optionally removes empty grid cells → reduces memory usage

Code

``` r

# Locate the example NetCDF file
ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")

# Read the NetCDF file with default settings
ncdata <- read_netcdf(ncfile)
```

#### Understanding the Output Structure

The
[`read_netcdf()`](https://github.com/Deltares-research/weathergenr/reference/read_netcdf.md)
function returns a list with five components that organize different
aspects of the climate data:

Code

``` r

# Display the top-level structure
names(ncdata)
```

    [1] "data"       "grid"       "date"       "dimensions" "attributes"

##### 1. Time Series Data

The `data` element contains a list of tibbles, one for each grid cell.
Each tibble represents a daily time series with columns for each
meteorological variable:

Code

``` r

# Examine the first grid cell's time series
ncdata$data[[1]]
```

    # A tibble: 7,671 × 7
       press_msl   kin temp_min temp_max  temp  kout precip
           <dbl> <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>
     1     1008.  174.     24.4     28.1  25.9  410.  9.12
     2     1008.  162.     23.7     28.1  24.9  410. 15.6
     3     1008.  201.     23.6     27.9  25.6  411.  4.84
     4     1008.  177.     25.0     28.3  26.1  411.  5.20
     5     1007.  218.     24.6     28.1  26.1  411.  1.54
     6     1007.  199.     25.0     28.2  26.3  411.  8.35
     7     1008.  195.     24.3     28.3  25.9  412. 14.4
     8     1009.  158.     24.6     26.5  25.4  412. 11.4
     9     1010.  143.     24.1     27.1  25.1  413.  9.15
    10     1009.  220.     23.9     28.6  26.1  413.  0.810
    # ℹ 7,661 more rows

This structure allows efficient access to complete time series for
individual locations while maintaining a compact memory footprint.

##### 2. Spatial Grid Information

The `grid` element provides spatial metadata for each grid cell:

Code

``` r

# View grid metadata
ncdata$grid
```

    # A tibble: 6 × 5
         id  xind  yind         x         y
      <int> <int> <int> <dbl[1d]> <dbl[1d]>
    1     1     1     1      9.5       0.5
    2     2     2     1      9.75      0.5
    3     3     3     1     10         0.5
    4     4     1     2      9.5       0.25
    5     5     2     2      9.75      0.25
    6     6     3     2     10         0.25

Each row describes a grid cell with:

- `id`: Sequential identifier for the grid cell
- `xind`, `yind`: Array indices in the original NetCDF file
- `x`, `y`: Spatial coordinates (longitude/latitude or projected
  coordinates)

##### 3. Temporal Information

The `date` element contains the complete time axis as a vector of Date
objects:

Code

``` r

# Display the first and last few dates
head(ncdata$date)
```

    [1] "2000-01-01" "2000-01-02" "2000-01-03" "2000-01-04" "2000-01-05"
    [6] "2000-01-06"

Code

``` r

tail(ncdata$date)
```

    [1] "2020-12-26" "2020-12-27" "2020-12-28" "2020-12-29" "2020-12-30"
    [6] "2020-12-31"

Code

``` r

# Confirm the total number of days
length(ncdata$date)
```

    [1] 7671

##### 4. Dimension Vectors

The `dimensions` element stores the raw dimension vectors from the
NetCDF file:

Code

``` r

# View available dimensions
names(ncdata$dimensions)
```

    [1] "time"      "longitude" "latitude" 

##### 5. Metadata Attributes

The `attributes` element preserves variable attributes, spatial
reference information, and global metadata from the NetCDF file,
ensuring traceability and documentation.

#### Advanced Usage

The
[`read_netcdf()`](https://github.com/Deltares-research/weathergenr/reference/read_netcdf.md)
function supports several optional parameters for customized data
loading:

Code

``` r

# Example: Load specific variables and rename them
ncdata <- read_netcdf(
  nc.file = ncfile,
  variables = c("precip", "temp"),
  var_rename = c(precip = "precipitation", temp = "temperature")
)

# Example: Handle 365-day calendar data (no leap days)
ncdata <- read_netcdf(
  nc.file = ncfile,
  leap.days = FALSE
)

# Example: Round values to reduce precision and memory
ncdata <- read_netcdf(
  nc.file = ncfile,
  signif.digits = 4,
  omit.empty = TRUE
)
```

## ⚙️ Generating Stochastic Weather Realizations

### Overview of the Weather Generation Method

The
[`generate_weather_series()`](https://github.com/Deltares-research/weathergenr/reference/generate_weather_series.md)
function implements a sophisticated multi-scale weather generation
framework that combines several statistical methods:

#### Statistical Components

1.  **Wavelet Autoregressive Modeling (WARM)**: Captures low-frequency
    (annual to multi-year) climate variability by decomposing annual
    precipitation signals into wavelet components → simulating their
    temporal evolution using ARIMA models.

2.  **Markov Chain Model**: Controls the persistence of wet and dry
    spells at the daily scale by modeling transitions between three
    precipitation states (dry, wet, extreme) with monthly-varying
    transition probabilities.

3.  **K-Nearest Neighbors (KNN) Resampling**: Selects historical
    analogue years based on similarity in annual precipitation patterns
    → resamples daily weather sequences from these analogues while
    preserving spatial and temporal correlations.

This coupled approach ensures that synthetic weather series preserve
both:

- **Temporal structure**: Day-to-day persistence, seasonal cycles, and
  interannual variability
- **Spatial coherence**: Correlations between grid cells and
  co-occurrence patterns of weather events
- **Statistical properties**: Means, variances, extremes, and
  distributional characteristics

The method is described in detail by [Steinschneider and Brown
(2013)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/wrcr.20528).

### Setting Up the Simulation

First, define the output directory and specify which meteorological
variables to simulate:

Code

``` r

# Create a temporary directory for outputs
output_path <- tempdir()

# Define variables to simulate
variables <- c("precip", "temp", "temp_min", "temp_max")

# Labels for plotting (optional)
variable_labels <- c("Precipitation", "Temperature (mean)", 
                     "Temperature (min)", "Temperature (max)")

# Number of stochastic realizations to generate
realization_num <- 3
```

#### Understanding Simulation Parameters

Key parameters control different aspects of the weather generation:

- **Temporal scope**: `sim.year.num` and `sim.year.start` define the
  simulation period
- **Realization count**: `realization.num` specifies how many
  independent synthetic series to create
- **WARM parameters**: Control annual-scale variability simulation
  - `warm.variable`: The primary driver variable (typically
    precipitation)
  - `warm.signif.level`: Significance threshold for retaining wavelet
    components (0.90 = 90% confidence)
  - `warm.sample.num`: Number of candidate annual traces to generate
    before filtering
- **KNN parameters**: Control daily-scale resampling
  - `knn.sample.num`: Number of historical analogue years to consider
- **Markov Chain parameters**: Define precipitation state thresholds
  - `mc.wet.quantile`: Threshold for wet day classification (0.3 = 30th
    percentile)
  - `mc.extreme.quantile`: Threshold for extreme precipitation (0.8 =
    80th percentile)

### Running the Weather Generator

Code

``` r

# Set a seed for reproducibility
seed <- 123

# Generate synthetic weather series
stochastic_weather <- generate_weather_series(
  weather.data = ncdata$data,
  weather.grid = ncdata$grid,
  weather.date = ncdata$date,
  variables = variables,
  sim.year.num = 20,
  sim.year.start = 2020,
  month.start = 1,
  realization.num = realization_num,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.num = 10000,
  knn.sample.num = 100,
  mc.wet.quantile = 0.2,
  mc.extreme.quantile = 0.8,
  output.path = output_path,
  seed = seed,
  compute.parallel = FALSE,
  num.cores = NULL
)
```

#### Understanding the Output

The function returns a list containing:

- `resampled`: The historical dates selected as analogues for each
  simulation day
- `dates`: The simulated daily time axis (always 365-day calendar, Feb
  29 excluded)

Code

``` r

# View the structure of the output
names(stochastic_weather)
```

    [1] "resampled" "dates"    

Code

``` r

# Check the simulated date range
head(stochastic_weather$dates)
```

    [1] "2020-01-01" "2020-01-02" "2020-01-03" "2020-01-04" "2020-01-05"
    [6] "2020-01-06"

Code

``` r

tail(stochastic_weather$dates)
```

    [1] "2039-12-26" "2039-12-27" "2039-12-28" "2039-12-29" "2039-12-30"
    [6] "2039-12-31"

Code

``` r

# Examine resampled dates for the first realization
head(stochastic_weather$resampled[[1]])
```

    [1] "2015-01-01" "2011-01-03" "2010-01-04" "2010-01-04" "2010-01-05"
    [6] "2010-01-06"

#### Interpreting Resampled Dates

The resampled dates show which historical days were selected as
analogues. This information reveals:

- Whether the generator is drawing from wet or dry historical periods
- The diversity of historical conditions being sampled
- Temporal clustering patterns in the resampling scheme

Code

``` r

# Display resampled dates for all realizations
stochastic_weather$resampled
```

    # A tibble: 7,300 × 3
       rlz_1      rlz_2      rlz_3
       <date[1d]> <date[1d]> <date[1d]>
     1 2015-01-01 2005-01-01 2013-01-01
     2 2011-01-03 2005-01-02 2013-01-05
     3 2010-01-04 2010-01-03 2013-01-06
     4 2010-01-04 2010-01-04 2013-01-02
     5 2010-01-05 2015-01-03 2004-01-05
     6 2010-01-06 2005-01-08 2004-01-07
     7 2010-01-10 2005-01-09 2004-01-08
     8 2010-01-11 2010-01-09 2004-01-25
     9 2010-01-12 2010-01-10 2013-01-13
    10 2010-01-13 2015-01-11 2004-01-08
    # ℹ 7,290 more rows

### Extracting Synthetic Weather Data

To access the actual weather values (not just the resampled dates), you
need to extract the corresponding data from the historical record:

Code

``` r

# Example: Extract synthetic weather for all realizations
synthetic_data <- list()

for (i in 1:realization_num) {
  # Get the day order for this realization
  day_order <- match(stochastic_weather$resampled[[i]], ncdata$date)
  
  # Extract corresponding weather data for each grid cell
  synthetic_data[[i]] <- lapply(ncdata$data, function(cell_data) {
    cell_data[day_order, variables] %>%
      mutate(date = stochastic_weather$dates, .before = 1)
  })
}
```

## 🔍 Evaluating Synthetic Weather Series

### Purpose of Evaluation

Validating the weather generator’s performance is crucial to ensure that
synthetic weather series are suitable for impact modeling. The
[`evaluate_weather_generator()`](https://github.com/Deltares-research/weathergenr/reference/evaluate_weather_generator.md)
function provides comprehensive diagnostic comparisons between simulated
and observed weather across multiple statistical dimensions.

### Preparing Data for Evaluation

The evaluation function requires:

1.  Observed weather data in list format (one element per grid cell)
2.  Simulated weather data structured similarly
3.  Matching date columns in both datasets

Code

``` r

# Create day order indices for each realization
day_order <- sapply(
  1:realization_num,
  function(n) match(stochastic_weather$resampled[[n]], ncdata$date)
)

# Extract synthetic weather samples
rlz_sample <- list()
for (n in 1:realization_num) {
  rlz_sample[[n]] <- lapply(ncdata$data[ncdata$grid$id], function(x) {
    x[day_order[, n], ] %>%
      select(precip, temp, temp_min, temp_max) %>%
      mutate(date = stochastic_weather$dates, .before = 1)
  })
}

# Prepare observed data
obs_sample <- lapply(ncdata$data[ncdata$grid$id], function(x) {
  x %>%
    select(precip, temp, temp_min, temp_max) %>%
    mutate(date = ncdata$date, .before = 1)
})
```

### Running the Evaluation

Code

``` r

# Evaluate weather generator performance
evaluation_plots <- evaluate_weather_generator(
  daily.sim = rlz_sample,
  daily.obs = obs_sample,
  output.path = NULL,
  save.plots = FALSE,
  variables = variables,
  variable.labels = variable_labels,
  realization.num = realization_num,
  wet.quantile = 0.3,
  extreme.quantile = 0.8
)
```

### Understanding Evaluation Metrics

The evaluation function computes and visualizes several key aspects of
weather generator performance:

#### 1. Daily Statistics

Comparison of mean, standard deviation, and skewness for each variable
across all grid cells. This assesses whether the generator preserves the
central tendency, variability, and distributional shape of observed
weather.

#### 2. Wet and Dry Day Frequencies

For precipitation, the evaluation examines:

- Total number of wet days (above the wet threshold)
- Number of extreme precipitation days (above the extreme threshold)
- Dry day frequencies

These metrics reveal whether the generator correctly reproduces
precipitation occurrence patterns.

#### 3. Spell Length Distributions

Analysis of consecutive wet and dry periods:

- Mean wet spell length
- Mean dry spell length
- Distribution of spell durations

This assesses the generator’s ability to capture persistence in weather
patterns.

#### 4. Spatial Correlations

Inter-site correlation matrices for each variable, comparing observed
versus simulated spatial coherence. This is critical for applications
requiring realistic spatial patterns (e.g., regional hydrology).

### Viewing Evaluation Results

Code

``` r

# Display all diagnostic plots
evaluation_plots
```

    $daily_mean

![](getting_started_files/figure-html/view_evaluation-1.png)


    $daily_sd

![](getting_started_files/figure-html/view_evaluation-2.png)


    $spell_length

![](getting_started_files/figure-html/view_evaluation-3.png)


    $wetdry_days_count

![](getting_started_files/figure-html/view_evaluation-4.png)


    $crossgrid

![](getting_started_files/figure-html/view_evaluation-5.png)


    $intergrid

![](getting_started_files/figure-html/view_evaluation-6.png)


    $precip_cond_cor

![](getting_started_files/figure-html/view_evaluation-7.png)


    $annual_pattern_precip

![](getting_started_files/figure-html/view_evaluation-8.png)


    $annual_pattern_temp

![](getting_started_files/figure-html/view_evaluation-9.png)


    $annual_pattern_temp_min

![](getting_started_files/figure-html/view_evaluation-10.png)


    $annual_pattern_temp_max

![](getting_started_files/figure-html/view_evaluation-11.png)


    $monthly_cycle

![](getting_started_files/figure-html/view_evaluation-12.png)


    $annual_precip

![](getting_started_files/figure-html/view_evaluation-13.png)


    attr(,"class")
    [1] "weather_assessment" "list"
    attr(,"fit_summary")
    # A tibble: 3 × 17
      rlz    rank overall_score mae_mean_precip mae_mean_temp mae_mean_temp_max
      <chr> <int>         <dbl>           <dbl>         <dbl>             <dbl>
    1 3         1         0.215            2.60          17.0              19.2
    2 1         2         0.498            2.78          17.0              19.2
    3 2         3         0.695            2.72          17.0              19.2
    # ℹ 11 more variables: mae_mean_temp_min <dbl>, mae_sd_precip <dbl>,
    #   mae_sd_temp <dbl>, mae_sd_temp_max <dbl>, mae_sd_temp_min <dbl>,
    #   mae_days_Dry <dbl>, mae_days_Wet <dbl>, mae_spell_Dry <dbl>,
    #   mae_spell_Wet <dbl>, mae_cor_crossgrid <dbl>, mae_cor_intervariable <dbl>
    attr(,"metadata")
    attr(,"metadata")$n_grids
    [1] 6

    attr(,"metadata")$realization.num
    [1] 3

    attr(,"metadata")$variables
    [1] "precip"   "temp"     "temp_min" "temp_max"

    attr(,"metadata")$assessment.date
    [1] "2026-01-11"

#### Interpreting the Plots

Good weather generator performance is indicated by:

- **Scatter plots near the 1:1 line**: Simulated statistics closely
  match observed values
- **Overlapping distributions**: Similar shapes and ranges between
  observed and simulated
- **Consistent correlation patterns**: Spatial correlations preserved in
  synthetic data
- **Balanced residuals**: No systematic bias in any direction

Discrepancies may indicate:

- Parameter tuning needed → adjust WARM or Markov chain thresholds
- Insufficient historical analogue pool → increase KNN sample size
- Structural limitations for specific climate regimes

## Other Features

The weathergenr package supports additional capabilities not covered in
this basic tutorial:

### Climate Scenario Analysis

Use
[`apply_climate_perturbations()`](https://github.com/Deltares-research/weathergenr/reference/apply_climate_perturbations.md)
to modify synthetic weather based on climate change scenarios:

Code

``` r

# Example: Apply temperature increase and precipitation changes
perturbed_weather <- apply_climate_perturbations(
  weather.data = synthetic_data[[1]],
  temperature.change = 2.0,  # 2°C increase
  precipitation.change = 1.1,  # 10% increase
  month.start = 1
)
```

### Wavelet Spectral Analysis

Explore low-frequency climate variability using wavelet transforms:

Code

``` r

# Analyze spectral characteristics of annual precipitation
wavelet_results <- wavelet_spectral_analysis(
  annual.series = annual_precip,
  signif.level = 0.90,
  mode = "fast"
)
```

### Parallel Computing

For large spatial domains, enable parallel processing to accelerate
computation:

Code

``` r

# Generate weather using multiple cores
stochastic_weather <- generate_weather_series(
  # ... other parameters ...
  compute.parallel = TRUE,
  num.cores = 4
)
```

### Saving and Exporting Results

Results can be written back to NetCDF format for use in other tools:

Code

``` r

# Example: Write synthetic weather to NetCDF
write_netcdf(
  weather.data = synthetic_data[[1]],
  weather.grid = ncdata$grid,
  weather.date = stochastic_weather$dates,
  variables = variables,
  output.file = file.path(output_path, "synthetic_weather.nc")
)
```

## Further Reading

- **Methodological details**: See Steinschneider and Brown (2013) for
  the theoretical foundation
- **Package documentation**: Use `?function_name` to access detailed
  help for any function

------------------------------------------------------------------------

**Session Information**

    R version 4.5.2 (2025-10-31 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26200)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=English_United States.utf8
    [2] LC_CTYPE=English_United States.utf8
    [3] LC_MONETARY=English_United States.utf8
    [4] LC_NUMERIC=C
    [5] LC_TIME=English_United States.utf8

    time zone: Europe/Amsterdam
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] ggplot2_4.0.1     dplyr_1.1.4       weathergenr_0.9.0

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6       xfun_0.55          lattice_0.22-7     quadprog_1.5-8
     [5] vctrs_0.6.5        tools_4.5.2        generics_0.1.4     curl_7.0.0
     [9] parallel_4.5.2     tibble_3.3.0       proxy_0.4-28       xts_0.14.1
    [13] pkgconfig_2.0.3    Matrix_1.7-4       RColorBrewer_1.1-3 S7_0.2.1
    [17] lifecycle_1.0.4    compiler_4.5.2     farver_2.1.2       textshaping_1.0.4
    [21] codetools_0.2-20   ncdf4_1.24         htmltools_0.5.9    class_7.3-23
    [25] yaml_2.3.12        pillar_1.11.1      tidyr_1.3.2        MASS_7.3-65
    [29] fitdistrplus_1.2-4 iterators_1.0.14   foreach_1.5.2      nlme_3.1-168
    [33] fracdiff_1.5-3     tidyselect_1.2.1   digest_0.6.39      purrr_1.2.0
    [37] labeling_0.4.3     tseries_0.10-58    splines_4.5.2      fastmap_1.2.0
    [41] grid_4.5.2         colorspace_2.1-2   cli_3.6.5          logger_0.4.1
    [45] magrittr_2.0.4     patchwork_1.3.2    utf8_1.2.6         survival_3.8-3
    [49] e1071_1.7-17       withr_3.0.2        scales_1.4.0       forecast_8.24.0
    [53] TTR_0.24.4         rmarkdown_2.30     quantmod_0.4.28    otel_0.2.0
    [57] nnet_7.3-20        timeDate_4051.111  gridExtra_2.3      ragg_1.5.0
    [61] zoo_1.8-15         urca_1.3-4         evaluate_1.0.5     knitr_1.51
    [65] lmtest_0.9-40      viridisLite_0.4.2  rlang_1.1.6        Rcpp_1.1.0
    [69] isoband_0.3.0      glue_1.8.0         rstudioapi_0.17.1  jsonlite_2.0.0
    [73] R6_2.6.1           systemfonts_1.3.1 
