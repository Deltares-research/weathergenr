
library(dplyr)
library(weathergenr)

rdata_path <- testthat::test_path("data", "testdata_generateWS.RData")
load(rdata_path)

realization_num <- rlz_num
stochastic_weather <- res
variables <- c("precip", "temp", "temp_min", "temp_max")
variable_labels <- variables


day_order <- sapply(
  1:realization_num,
  function(n) match(stochastic_weather$resampled[[n]], ncdata$date)
)

rlz_sample <- list()
for (n in 1:realization_num) {
  rlz_sample[[n]] <- lapply(ncdata$data[ncdata$grid$id], function(x) {
    x[day_order[, n], ] %>%
      select(precip, temp, temp_min, temp_max) %>%
      mutate(date = stochastic_weather$dates, .before = 1)
  })
}

obs_sample <- lapply(ncdata$data[ncdata$grid$id], function(x) {
  x %>%
    select(precip, temp, temp_min, temp_max) %>%
    dplyr::mutate(date = ncdata$date, .before = 1)
})

daily.sim <- rlz_sample
daily.obs <- obs_sample
save.plots <- TRUE
variables <- variables
variable.labels <- variable_labels
variable.units <- NULL
realization.num <- realization_num
wet.quantile <- 0.3
extreme.quantile <- 0.8
output.path = "C:/TEMP/EVAL/"


timez <- Sys.time()

out <- evaluateWegen(
  output.path = output.path,
  daily.sim = rlz_sample,
  daily.obs = obs_sample,
  save.plots = save.plots,
  variables = variables,
  variable.labels = variable_labels,
  variable.units = NULL,
  realization.num = realization_num,
  wet.quantile = 0.3,
  extreme.quantile = 0.8
)
Sys.time() - timez
