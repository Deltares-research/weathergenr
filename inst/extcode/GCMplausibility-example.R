
# Load GCM projections
gcm_data <- readr::read_csv("./inst/extdata/annual_change_scalar_stats_summary_mean.csv")
gcm_data <- gcm_data %>% filter(horizon %in% "near") %>% na.omit()
gcm_data <- gcm_data %>% select(scenario, model, x = prcp, y = tavg)

# Load stress test results
str_data <- readr::read_csv("./inst/extdata/Qstats.csv")

# Run function
gcm_plausibility <- GCMplausiblity(str.data = str_data,
                                   gcm.data = gcm_data,
                                   clevel.list = c(0.50, 0.95),
                                   gcm.scenario.list = c("ssp126", "ssp245", "ssp370", "ssp585"),
                                   metric.list = c("mean","max","min", "q95", "Q7day_min"),
                                   metric.labs = NULL,
                                   location.list = c("Q_1"))

# Save results to csv
readr::write_csv(gcm_plausibility, file = "./temp/gcm_plausibility.csv")


