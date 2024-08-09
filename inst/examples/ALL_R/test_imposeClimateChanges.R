

testthat::test_that("check if the output dimensions are correct", {

  # Calculated perturbed statistics
  output <- imposeClimateChanges(climate.data = stochastic_rlz,
                                 climate.grid = ncdata$grid,
                                 sim.dates = stochastic_weather$dates,
                                 change.factor.precip.mean = delta_precip_mean,
                                 change.factor.precip.variance = delta_precip_variance,
                                 change.factor.temp.mean = delta_temp_mean,
                                 transient.temp.change = TRUE,
                                 transient.precip.change = TRUE,
                                 calculate.pet = TRUE,
                                 compute.parallel = TRUE,
                                 num.cores = NULL,
                                 fit.method = "mme"
  )

  # Compare against baseline

})

