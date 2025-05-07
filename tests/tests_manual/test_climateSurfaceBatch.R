
#devtools::install_github("deltares/weathergenr", ref = "Dev")

################################################################################

source("./tests/tests_manual/test_loadData.R")

########################################

metrics_to_plot <- c("mean", "min", "max", "q95", "Q7day_max", "Q7day_min")
locations_to_plot <- c("Q_1", "Q_5001", "Q_4001")

climateSurfaceBatch(
  str.data = str_data_rel,
  gcm.data = gcm_data,
  save.dir = sdir,
  variable.x.breaks = seq(-30,30,10),
  variable.y.breaks = seq(0,3,1),
  gcm.legend = TRUE,
  gcm.future.period = "near",
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = c(0.5, 0.75, 0.95),
  metrics = metrics_to_plot,
  metric.labels = paste0("Metric: ", metrics_to_plot),
  locations = locations_to_plot,
  location.labels = paste0("Location: ", locations_to_plot),
  metric.direction = c(1,1,0,1,0,1,1,0),
  relative.results = FALSE,
  text.scale = 0.6)


