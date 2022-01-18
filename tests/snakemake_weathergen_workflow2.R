
yaml2 <- yaml::read_yaml("stresstest_config.yml")

# General attributes
output_path <- yaml2$general$output.path
ncfile_prefix <- yaml2$general$nc.file.prefix
nmax <- yaml2$general$historical_realizations_num

# Temperature change attributes
delta_temp_mean_min <-  yaml2$temp$mean$min
delta_temp_mean_max <-  yaml2$temp$mean$max
temp_step_num <- yaml2$temp$step_num
temp_change_type <- yaml2$temp$change_type


# Precip change attributes
delta_precip_mean_min <- yaml2$precip$mean$min
delta_precip_mean_max <- yaml2$precip$mean$max
delta_precip_variance_min <- yaml2$precip$variance$min
delta_precip_variance_max <- yaml2$precip$variance$max
precip_step_num <- yaml2$precip$step_num
precip_change_type <- yaml2$precip$change_type

# Prepare stress test matrix

temp_mean_steps <- sapply(1:12, function(m)
         seq(delta_temp_mean_min[m], delta_temp_mean_max[m],
             length.out = temp_step_num))

precip_mean_steps <- sapply(1:12, function(m)
         seq(delta_precip_mean_min[m], delta_precip_mean_max[m],
             length.out = precip_step_num))

precip_variance_steps <- sapply(1:12, function(m)
         seq(delta_precip_variance_min[m], delta_precip_variance_max[m],
             length.out = precip_step_num))

strtest_matrix <- tidyr::expand_grid(precip_ind = 1:precip_step_num, temp_ind = 1:temp_step_num)
 write.csv(strtest_matrix,
   paste0(output_path, "strtest_matrix.csv"), row.names = FALSE)


# Stress test delta factors for each variable/climate statistic
strtest_matrix_precip_mean <- precip_mean_steps[strtest_matrix$precip_ind, ]
 write.csv(strtest_matrix_precip_mean,
   paste0(output_path, "strtest_matrix_precip_mean.csv"), row.names = FALSE)

strtest_matrix_precip_variance <- precip_variance_steps[strtest_matrix$precip_ind, ]
 write.csv(strtest_matrix_precip_variance,
   paste0(output_path, "strtest_matrix_precip_variance.csv"), row.names = FALSE)

strtest_matrix_temp_mean <- temp_mean_steps[strtest_matrix$temp_ind, ]
 write.csv(strtest_matrix_temp_mean,
   paste0(output_path, "strtest_matrix_temp_mean.csv"), row.names = FALSE)

# Total climate changes
smax <- nrow(strtest_matrix)

hist_rlz_vector <- list.files("C:/testrun2/historical", pattern = ".nc", full.names = TRUE)

# Loop through weather realizations
for (n in 1:nmax) {

  rlz_input <- readNetcdf(hist_rlz_vector[n])

  for (s in 1:smax) {

    # Apply climate changes to climate data
    rlz_future <- imposeClimateChanges(
       climate.data = rlz_input$data,
       climate.grid = rlz_input$grid,
       sim.dates = rlz_input$date,
       change.factor.precip.mean = strtest_matrix_precip_mean[s,],
       change.factor.precip.variance = strtest_matrix_precip_variance[s,],
       change.factor.temp.mean = strtest_matrix_temp_mean[s,],
       change.type.temp = temp_change_type,
       change.type.precip = precip_change_type)

     # Save to netcdf file
     writeNetcdf(
       data = rlz_future,
       coord.grid = rlz_input$grid,
       output.path = output_path,
       origin.date =  rlz_input$date,
       calendar.type = "noleap",
       nc.template.file = hist_rlz_vector[n],
       nc.compression = 4,
       nc.spatial.ref = "spatial_ref",
       nc.file.prefix = yaml2$general$nc.file.prefix,
       nc.file.suffix = paste0(n,"_",s))
  }
}

