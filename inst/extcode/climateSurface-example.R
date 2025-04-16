
################################################################################

devtools::install_github("deltares/weathergenr", ref = "Dev")

library(weathergenr)
library(dplyr)
library(ggplot2)

# Load GCM projections
gcm_data <- readr::read_csv("./tests/testdata/annual_change_scalar_stats_summary_mean.csv")
gcm_data <- gcm_data %>% filter(horizon %in% "near")

# Load stress test results
str_data <- readr::read_csv("./tests/testdata/Qstats.csv") %>% filter(statistic == "min")
sdir <- "./temp/"

################################################################################

# Calculate relative changes
bindex <- which(str_data$tavg == 0 & str_data$prcp == 0)
str_data_rel <- str_data %>% mutate(Q_1 = Q_1/str_data[[bindex,"Q_1"]] * 100 - 100)

# Plot response surface
px <- climateSurface(
  str.data = str_data_rel,
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
  gcm.bvnorm.levels = c(0.5, 0.95),
  gcm.transparency = 0.75,
  gcm.legend = TRUE,
  variable.x.breaks = seq(-40, 40, 10),
  variable.y.breaks = seq(0,6,1),
  variable.z.breaks = NULL, #seq(-160, 260, 40),
  text.scale = 0.6,
  contour.num = 20,
  multi.panel = FALSE,
  panel.variable = NULL,
  panel.variable.levels = NULL)

#Save plot to file
ggsave(filename=paste0(sdir,"crsplot1.png"), height = 6.5, width = 7.5, dpi = 500)

