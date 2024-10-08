## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE, cache = TRUE)

## ----eval, eval = TRUE, cache=TRUE--------------------------------------------
weather_data <- ncdata$data

day_order <- sapply(1:realization_num,
      function(n) match(stochastic_weather$resampled[[n]], ncdata$date))

rlz_sample <- list()
for (n in 1:realization_num) {
  rlz_sample[[n]] <- lapply(weather_data[ncdata$grid$id], function(x) x[day_order[,n],] %>%
                              select(precip,  temp, temp_min, temp_max) %>%
                              mutate(date = stochastic_weather$dates, .before = 1))
}

obs_sample <- lapply(weather_data[ncdata$grid$id], function(x) x %>%
                       select(precip, temp, temp_min, temp_max) %>%
                       dplyr::mutate(date = ncdata$date, .before = 1))

out <- evaluateWegen(daily.sim = rlz_sample,
              daily.obs = obs_sample,
              output.path = paste0(output_path,"/plots"),
              variables = c("precip", "temp", "temp_min", "temp_max"),
              variable.labels = c("Precip.", "Temp. (avg)", "Temp. (min)", "Temp. (max)"),
              variable.units = NULL,
              realization.num = realization_num,
              wet.quantile = 0.2,
              extreme.quantile =0.8)


