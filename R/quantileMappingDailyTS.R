


quantileMappingDailyTS <- function(
  parms = NULL,
  base.ts = NULL,
  base.dist = "gamma", 
  change.operator = "multiply") {

  require(dplyr)
  require(readr)
  require(fitdistrplus)
  require(lubridate)
 
  # Number of years in the time-series
  ymax <-  max(year(base.ts$date)) - year(base.ts$date[1]) + 1
  
  # Empty lists and vectors to store distributions
  emp1 <- vector("list", length = 12)
  bdist <- list(fit = emp1, shape = emp1, scale = emp1, mean = emp1, var = emp1)
  emp2 <- matrix(0, nrow = 12, ncol = ymax)
  tdist <- list(fit = emp2, shape = emp2, scale = emp2, mean = emp2, var = emp2)
  
  # Identify months to be perturbed
  if(change.operator == "mult") {nochange_mean <- 1; nochange_var <- 1 
  } else {nochange_mean <- 0; nochange_var <- 1}
  
  # Estimate historical quantiles and map to the target distribution 
  base.ts$mon  <- month(base.ts$date)
  base.ts$year <- year(base.ts$date) - min(year(base.ts$date)) + 1
  
  #Keep track of changing months
  pmon <- which((parms$mean != nochange_mean) | (parms$var != nochange_var))  
  pind <- which(base.ts$mon %in% pmon)
  
  if (length(pind) == 0) {
    return(base.ts$value)
    
  } else {
    
    new_ts <- base.ts[pind, ]
    new_ts_nonzero <- new_ts[new_ts$value > 0,]
    
    # Fit base distribution to shifting months
    bdist[["fit"]][pmon] <- lapply(pmon, function(x) quiet(fitdist(new_ts_nonzero %>% 
                                                                     filter(month(date) == x) %>% pull(value), distr = base.dist, method = "mle")))
    
    bdist[["shape"]][pmon] <- sapply(pmon, function(x) bdist[["fit"]][[x]]$estimate[[1]])
    bdist[["scale"]][pmon] <- sapply(pmon, function(x) 1/(bdist[["fit"]][[x]]$estimate[[2]]))
    bdist[["mean"]][pmon]  <- sapply(pmon, function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]])
    bdist[["var"]][pmon]   <- sapply(pmon, function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]]^2)
    
    if (change.operator == "multiply") {
      
      # Vector of means and variances (row=years, column=months)
      parms_yearly <- list()
      parms_yearly[["mean"]] <- sapply(1:12, 
                                       function(m) seq(1, parms$mean[m], length.out = ymax))
      parms_yearly[["var"]]  <- sapply(1:12, 
                                       function(m) seq(1, parms$var[m], length.out = ymax))
      
      # Define the mean, variance of the target distribution
      tdist[["mean"]][pmon,] <- sapply(1:ymax, function(y) sapply(pmon, 
                                                                  function(m) bdist[["mean"]][[m]] * parms_yearly[["mean"]][y, m]))
      
      tdist[["var"]][pmon,] <- sapply(1:ymax, function(y) sapply(pmon, 
                                                                 function(m) bdist[["var"]][[m]] * parms_yearly[["var"]][y, m]))
      
    } else if (change.operator == "add") {
      
      # Vector of means and variances (row=years, column=months)
      parms_yearly <- list()
      parms_yearly[["mean"]] <- sapply(1:12, 
                                       function(m) seq(0, parms$mean[m], length.out = ymax))
      parms_yearly[["var"]]  <- sapply(1:12, 
                                       function(m) seq(1, parms$var[m], length.out = ymax))
      
      # Define the mean, variance of the target distribution
      tdist[["mean"]][pmon,] <- sapply(1:ymax, function(y) sapply(pmon, 
                                                                  function(m) bdist[["mean"]][[m]] + parms_yearly[["mean"]][y, m]))
      
      tdist[["var"]][pmon,] <- sapply(1:ymax, function(y) sapply(pmon, 
                                                                 function(m) bdist[["var"]][[m]] * parms_yearly[["var"]][y, m]))
      
    } else {
      stop("change operator can be either 'multiply' or 'add'")
    }
    
    # Define the shape and scale parameters of the target distribution
    tdist[["scale"]][pmon,] <- sapply(1:ymax, function(y) sapply(pmon, 
                                                                 function(m) tdist[["var"]][m,y] / tdist[["mean"]][m,y]))
    
    tdist[["shape"]][pmon,] <- sapply(1:ymax, function(y) sapply(pmon, 
                                                                 function(m) tdist[["mean"]][m,y]^2 / tdist[["var"]][m,y]))
    
    
    qtile <- pgamma(new_ts$value, shape = unlist(bdist[["shape"]][new_ts$mon]), 
                    scale = unlist(bdist[["scale"]][new_ts$mon]))
    
    shape_vec <- sapply(1:nrow(new_ts), function(x) tdist[["shape"]][new_ts$mon[x],new_ts$year[x]])
    scale_vec <- sapply(1:nrow(new_ts), function(x) tdist[["scale"]][new_ts$mon[x],new_ts$year[x]])
    
    # Replace the original matrix
    base.ts$value[pind] <- qgamma(qtile, shape = shape_vec, scale = scale_vec)
    
    return(base.ts$value)

  }
  
  
 
} 