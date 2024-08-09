
## ----setup, eval=TRUE, echo = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(weathergenr)
library(dplyr)
library(tidyr)
library(ggplot2)


## ----load-data---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
ncdata <- readNetcdf(ncfile)

resampled_dates <- readr::read_csv("../tests/testdata/resampled_dates.csv")
simulated_dates <- readr::read_csv("../tests/testdata/sim_dates.csv")
realization.num <- ncol(resampled_dates)


## ----parameters--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Maximum number of grid cells to consider for performance assessment
evaluate.grid.num <- 20


## ----prepare-files-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grids <- ncdata$grid$id
sampleGrids <- sample(grids, size = min(evaluate.grid.num, length(grids)))

daily_obs <- lapply(ncdata$data[sampleGrids], function(x)
        dplyr::mutate(x, date = ncdata$date, .before = 1))

daily_sim <- vector(mode = "list", length = realization.num)

for (n in 1:realization.num) {

  # Resample order
  day_order <- match(resampled_dates[[n]], ncdata$date)

  # Obtain stochastic series by re-ordering historical data
  daily_sim[[n]] <- lapply(ncdata$data[sampleGrids], function(x)
    x[day_order,] %>% mutate(date = simulated_dates[[1]], .before = 1))
}


## ----run-evaluate------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plots <- evaluateWegen(daily.sim = daily_sim,
              daily.obs = daily_obs,
              output.path = tempdir(),
              variables = c("precip", "temp", "temp_min", "temp_max"),
              variable.labels = c("precip", "temp", "temp_min", "temp_max"),
              variable.units = NULL,
              realization.num = realization.num,
              wet.quantile = 0.2,
              extreme.quantile = 0.8,
              show.title = TRUE)



## ----run-visualize, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------
print(plots[1])

