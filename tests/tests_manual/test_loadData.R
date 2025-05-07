
# Packages required
library(weathergenr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load GCM projections
gcm_data <- readr::read_csv("./tests/testdata/annual_change_scalar_stats_summary_mean.csv") %>%
  rename(prcp = precip, tavg = temp)
gcm_data <- gcm_data %>% filter(horizon %in% "near")

# Load stress test results
str_data <- readr::read_csv("./tests/testdata/Qstats.csv")
sdir <- "./tests/temp/"
file.remove(list.files(path = sdir, full.names = TRUE))

# Calculate relative changes
str_data_rel <- str_data %>%
  pivot_longer(cols = starts_with("Q_"), names_to = "location", values_to = "value") %>%
  group_by(statistic, location) %>%
  mutate(value = value/.data[["value"]][which(.data$tavg == 0 & .data$prcp == 0)] * 100 - 100) %>%
  pivot_wider(names_from = "location", values_from = "value")
