utils::globalVariables(c(
  "scenario", "prcp", "tavg", "statistic",
  "z", "y", "clevel",
  "Location", "Metric", "Baseline", "colorRampPalette",".",
  'level_mid', 'horizon',
  'value', 'minval', 'maxval'

))

#Workaround for rlang warning
year <- mon <- day <- precip <- sd <- variable <- value <- id_variable <- 0
Wet <- Dry <- wet <- dry <- id_variable1 <- id_variable2 <- 0
rlz <- id1 <- id2 <- variable1 <- variable2 <- type <- Observed <- Stochastic <- 0
wet_th <- Simulated <- 0
