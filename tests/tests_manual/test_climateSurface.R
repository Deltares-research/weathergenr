
#devtools::install_github("deltares/weathergenr", ref = "Dev")

################################################################################

source("./tests/tests_manual/test_loadData.R")

########################################

metric <- "mean"

# Plot response surface only
px <- (climateSurface(
  str.data = str_data_rel %>% filter(statistic == metric),
  gcm.data = NULL,
  variable.x = "prcp",
  variable.y = "tavg",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = bquote(bold('Change in mean flow (%)')),
  variable.x.label = expression(Delta~"Precipitation"),
  variable.y.label = expression(Delta~"Temperature"),
  failure.direction = 1,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = NULL,
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL, #seq(-160, 260, 40),
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)) %>%
ggsave(filename=paste0(sdir,"crsplot1.png"), height = 6.5, width = 7.5, dpi = 500)


# Plot response surface only + Default z-range
px <- (climateSurface(
  str.data = str_data_rel %>% filter(statistic == metric),
  gcm.data = NULL,
  variable.x = "prcp",
  variable.y = "tavg",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = bquote(bold('Change in mean flow (%)')),
  variable.x.label = expression(Delta~"Precipitation"),
  variable.y.label = expression(Delta~"Temperature"),
  failure.direction = 1,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = NULL,
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL,
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)) %>%
  ggsave(filename=paste0(sdir,"/crsplot2.png"), height = 6.5, width = 7.5, dpi = 500)



# Plot response surface only + Default z-range + change direction
px <- (climateSurface(
  str.data = str_data_rel %>% filter(statistic == metric),
  gcm.data = NULL,
  variable.x = "prcp",
  variable.y = "tavg",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = bquote(bold('Change in mean flow (%)')),
  variable.x.label = expression(Delta~"Precipitation"),
  variable.y.label = expression(Delta~"Temperature"),
  failure.direction = 0,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = NULL,
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL,
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)) %>%
  ggsave(filename=paste0(sdir,"/crsplot3.png"), height = 6.5, width = 7.5, dpi = 500)


# Plot response surface + GCM Dots
px <- (climateSurface(
  str.data = str_data_rel %>% filter(statistic == metric),
  gcm.data = gcm_data,
  variable.x = "prcp",
  variable.y = "tavg",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = bquote(bold('Change in mean flow (%)')),
  variable.x.label = expression(Delta~"Precipitation"),
  variable.y.label = expression(Delta~"Temperature"),
  failure.direction = 1,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = NULL,
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL,
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)) %>%
  ggsave(filename=paste0(sdir,"/crsplot4.png"), height = 6.5, width = 7.5, dpi = 500)


# Plot response surface + GCM Dots + no legend
px <- (climateSurface(
  str.data = str_data_rel %>% filter(statistic == metric),
  gcm.data = gcm_data,
  variable.x = "prcp",
  variable.y = "tavg",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = bquote(bold('Change in mean flow (%)')),
  variable.x.label = expression(Delta~"Precipitation"),
  variable.y.label = expression(Delta~"Temperature"),
  failure.direction = 1,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = NULL,
  gcm.transparency = 0.75,
  gcm.legend = FALSE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL,
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)) %>%
  ggsave(filename=paste0(sdir,"/crsplot5.png"), height = 6.5, width = 7.5, dpi = 500)

# Plot response surface + GCM Dots + confidence interval
px <- (climateSurface(
  str.data = str_data_rel %>% filter(statistic == metric),
  gcm.data = gcm_data,
  variable.x = "prcp",
  variable.y = "tavg",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = bquote(bold('Change in mean flow (%)')),
  variable.x.label = expression(Delta~"Precipitation"),
  variable.y.label = expression(Delta~"Temperature"),
  failure.direction = 1,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = c(0.5, 0.75, 0.99),
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL,
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)) %>%
  ggsave(filename=paste0(sdir,"/crsplot6.png"), height = 6.5, width = 7.5, dpi = 500)


#############################################################################

# Climate Surface Panel
metrics_to_plot <- c("mean", "min", "max", "q95", "Q7day_max", "Q7day_min")

panel_row_num <- ceiling(length(metrics_to_plot)/2)



# Plot response surface + GCM Dots + confidence interval
px <- (climateSurface(
  str.data = str_data_rel,
  gcm.data = gcm_data,
  variable.x = "prcp",
  variable.y = "tavg",
  variable.z = "Q_1",
  threshold.z = NULL,
  plot.title = bquote(bold('Change in metrics (%)')),
  variable.x.label = expression(Delta~"Precipitation"),
  variable.y.label = expression(Delta~"Temperature"),
  failure.direction = 1,
  gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
  gcm.bvnorm.levels = c(0.5, 0.75, 0.99),
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL,
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = TRUE,
  panel.variable = "statistic",
  panel.variable.levels = metrics_to_plot)) %>%
  ggsave(filename=paste0(sdir,"/crsplot7.png"),
         height = 6.5 + 2.5 * (panel_row_num-1),
         width = 7.5, dpi = 500)


