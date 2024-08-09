
library(dplyr)
library(tidyr)
library(weathergenr)
library(ggplot2)


evaluate.grid.num <- 25

# Data
plots_path <- "C:/Users/taner/Workspace/test-weathergenr/data/vakhsh/"
ncfile <- paste0(plots_path, "extract_historical.nc")
sim_dates <- readr::read_csv(paste0(plots_path, "sim_dates.csv"))
resampled_dates <- readr::read_csv(paste0(plots_path, "resampled_dates.csv"))

ncdata <- weathergenr::readNetcdf(ncfile)
realization.num <- ncol(resampled_dates)
ngrids <- length(ncdata$data)
weather.date <- ncdata$date
weather.data <- ncdata$data
weather.grid <- ncdata$grid
sampleGrids <- sample(grids, size = min(evaluate.grid.num, ngrids))


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

evaluateWegen(daily.sim = rlz_sample,
                daily.obs = obs_sample,
                output.path = paste0(plots_path,"v3/"),
                variables = c("precip", "temp", "temp_min", "temp_max"),
                variable.labels = c("Precip.", "Temp. (avg)", "Temp. (min)", "Temp. (max)"),
                variable.units = NULL,
                realization.num = realization.num,
                wet.quantile = 0.3,
                extreme.quantile =0.8)



