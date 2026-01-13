
utils::globalVariables(c(
  "scenario", "prcp", "tavg", "statistic",
  "z", "y", "clevel",
  "Location", "Metric", "Baseline", "colorRampPalette", ".",
  "level_mid", "horizon",
  "value", "minval", "maxval",
  "min_val", "max_val", "lim", "minlim", "maxlim",
  "pct.change","assessment","original","adjusted","threshold",
  "idx","series","period","lo","hi","obs","sim_mu",
  "xmin","xmax","ymin","ymax","dateo"
))

# Workaround for rlang warning
wyear <- year <- mon <- day <- precip <- sd <- variable <- value <- id_variable <- 0
Wet <- Dry <- wet <- dry <- id_variable1 <- id_variable2 <- 0
rlz <- id1 <- id2 <- variable1 <- variable2 <- type <- Observed <- Stochastic <- 0
wet_th <- Simulated <- 0
x <- y <- z <- 0
metric <- baseline <- low <- high <- month <- 0
