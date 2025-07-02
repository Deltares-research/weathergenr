

source("./tests/tests_manual/test_loadData.R")

# Select the metrics to include in the calculation
metrics_list <- c("mean", "median", "q70")

# Set labels for the selected metrics
metrics_labels <- c('Mean Flow (%)', 'Median Flow (%)', 'Q70 Discharge (%)')

# Select discharge locations for the analysis
location_list <- c("Q_5001", "Q_5002")

# Run function
gcm_plausibility <- GCMplausiblity(str.data = str_data,
                                   gcm.data = gcm_data,
                                   gcm.scenario.list = c("ssp126", "ssp245", "ssp370", "ssp585"),
                                   clevel.list = c(0.50, 0.95),
                                   metric.list = metrics_list,
                                   metric.labs = metrics_labels,
                                   location.list = location_list)

#readr::write_csv(gcm_plausibility, file = paste0(sdir, "gcm_plausibility.csv"))
