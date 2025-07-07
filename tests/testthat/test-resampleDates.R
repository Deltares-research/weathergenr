
library(testthat)


# test_that("check resampleDates returns same results with same seed", {
#
#   rdata_path <- testthat::test_path("data", "testdata_resampleDates.RData")
#   e <- new.env()
#   load(rdata_path, envir = e)
#
#   seed <- 100
#
#   sim1 <- resampleDates(
#     PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
#     ANNUAL_PRCP = warm_variable,
#     PRCP = climate_d_aavg$precip,
#     TEMP = climate_d_aavg$temp,
#     TMAX = climate_d_aavg$temp_max,
#     TMIN = climate_d_aavg$temp_min,
#     START_YEAR_SIM = sim.year.start,
#     k1 = n,
#     ymax = sim.year.num,
#     dates.d = dates_d,
#     sim.dates.d = sim_dates_d,
#     knn.annual.sample.num = knn.sample.num,
#     dry.spell.change = dry.spell.change,
#     wet.spell.change = wet.spell.change,
#     YEAR_D = year_seq,
#     month.start = month.start,
#     wet.quantile = mc.wet.quantile,
#     extreme.quantile = mc.extreme.quantile,
#     seed = seed + n)
#
#   sim2 <- resampleDates(
#     PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
#     ANNUAL_PRCP = warm_variable,
#     PRCP = climate_d_aavg$precip,
#     TEMP = climate_d_aavg$temp,
#     TMAX = climate_d_aavg$temp_max,
#     TMIN = climate_d_aavg$temp_min,
#     START_YEAR_SIM = sim.year.start,
#     k1 = n,
#     ymax = sim.year.num,
#     dates.d = dates_d,
#     sim.dates.d = sim_dates_d,
#     knn.annual.sample.num = knn.sample.num,
#     dry.spell.change = dry.spell.change,
#     wet.spell.change = wet.spell.change,
#     YEAR_D = year_seq,
#     month.start = month.start,
#     wet.quantile = mc.wet.quantile,
#     extreme.quantile = mc.extreme.quantile,
#     seed = seed + n)
#
#   expect_identical(sim1, sim2)
# })
#
# test_that("check resampleDates returns different results with different seeds", {
#
#   load("./tests/testdata/testdata_resampleDates.RData")
#   seed <- 100
#
#   sim1 <- resampleDates(
#     PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
#     ANNUAL_PRCP = warm_variable,
#     PRCP = climate_d_aavg$precip,
#     TEMP = climate_d_aavg$temp,
#     TMAX = climate_d_aavg$temp_max,
#     TMIN = climate_d_aavg$temp_min,
#     START_YEAR_SIM = sim.year.start,
#     k1 = n,
#     ymax = sim.year.num,
#     dates.d = dates_d,
#     sim.dates.d = sim_dates_d,
#     knn.annual.sample.num = knn.sample.num,
#     dry.spell.change = dry.spell.change,
#     wet.spell.change = wet.spell.change,
#     YEAR_D = year_seq,
#     month.start = month.start,
#     wet.quantile = mc.wet.quantile,
#     extreme.quantile = mc.extreme.quantile,
#     seed = seed + n)
#
#   seed <- 200
#
#   sim2 <- resampleDates(
#     PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
#     ANNUAL_PRCP = warm_variable,
#     PRCP = climate_d_aavg$precip,
#     TEMP = climate_d_aavg$temp,
#     TMAX = climate_d_aavg$temp_max,
#     TMIN = climate_d_aavg$temp_min,
#     START_YEAR_SIM = sim.year.start,
#     k1 = n,
#     ymax = sim.year.num,
#     dates.d = dates_d,
#     sim.dates.d = sim_dates_d,
#     knn.annual.sample.num = knn.sample.num,
#     dry.spell.change = dry.spell.change,
#     wet.spell.change = wet.spell.change,
#     YEAR_D = year_seq,
#     month.start = month.start,
#     wet.quantile = mc.wet.quantile,
#     extreme.quantile = mc.extreme.quantile,
#     seed = seed + n)
#
#   expect_false(identical(sim1, sim2))
#
# })
#
# test_that("check resampleDates handles short records", {
#
#   load("./tests/testdata/testdata_resampleDates.RData")
#
#   seed <- 100
#   sp <- 1:(1096*2) #6years
#   yrnm <- floor(max(sp)/365)
#
#   sim1 <- resampleDates(
#     PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
#     ANNUAL_PRCP = warm_variable[1:yrnm],
#     PRCP = climate_d_aavg$precip[sp],
#     TEMP = climate_d_aavg$temp[sp],
#     TMAX = climate_d_aavg$temp_max[sp],
#     TMIN = climate_d_aavg$temp_min[sp],
#     START_YEAR_SIM = sim.year.start,
#     k1 = n,
#     ymax = sim.year.num,
#     dates.d = dates_d[sp,],
#     sim.dates.d = sim_dates_d,
#     knn.annual.sample.num = knn.sample.num,
#     dry.spell.change = dry.spell.change,
#     wet.spell.change = wet.spell.change,
#     YEAR_D = year_seq[sp],
#     month.start = month.start,
#     wet.quantile = mc.wet.quantile,
#     extreme.quantile = mc.extreme.quantile,
#     seed = seed + n)
#
# })
