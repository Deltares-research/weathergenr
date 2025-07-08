

library(weathergenr)
library(testthat)
rdata_path <- testthat::test_path("data", "testdata_resampleDates.RData")
load(rdata_path)
n=1

#Global seed
set.seed(999)

test_that("check results are generated", {

  cur_time <- Sys.time()

  sim1 <- resampleDates(
    PRCP_FINAL_ANNUAL_SIM = sim_annual_sub$sampled[, n],
    ANNUAL_PRCP = warm_variable,
    PRCP = climate_d_aavg$precip,
    TEMP = climate_d_aavg$temp,
    START_YEAR_SIM = sim.year.start,
    k1 = n,
    ymax = sim.year.num,
    dates.d = dates_d,
    sim.dates.d = sim_dates_d,
    knn.annual.sample.num = knn.sample.num,
    dry.spell.change = dry.spell.change,
    wet.spell.change = wet.spell.change,
    YEAR_D = year_seq,
    month.start = month.start,
    wet.quantile = mc.wet.quantile,
    extreme.quantile = mc.extreme.quantile,
    seed = 100)

  print(Sys.time() - cur_time)

  expect_true(length(sim1) == length(sim_dates_d$date))
  #expect_true(sum(as.numeric(sim1)) == 107283671) # with seed 100
})

