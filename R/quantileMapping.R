
#' Function for perturbing climate statistics using Quantile-Mapping
#'
#' @param par List of parameters to be passed.
#' @param operator A character string to specify the type of perturb change operator.
#' @param value To be completed...
#' @param date  To be completed...
#'
#' @return
#' @export
#' @importFrom fitdistrplus fitdist
#'
quantileMapping <- function(
  value = NULL,
  date = NULL,
  par = NULL,
  mon.ts = NULL,
  year.ts = NULL)

  {

  emp1 <- matrix(0, nrow = 12, ncol = max(year.ts))
  emp2 <- vector("list", 12)

  bdist <- list(fit = emp2, shape = emp2, scale = emp2, mean = emp2, var = emp2)
  tdist <- list(fit = NA, shape = emp1, scale = emp1, mean = emp1, var = emp1)

  # Identify months to be perturbed
  nochange_var = 1
  nochange_mean <- 1

  # Vector of means and variances (row=years, column=months)
  parms_yearly <- list()
  parms_yearly[["mean"]] <- sapply(1:12, function(m) seq(nochange_mean, par$mean[m], length.out = max(year.ts)))
  parms_yearly[["var"]]  <- sapply(1:12, function(m) seq(nochange_var, par$var[m], length.out = max(year.ts)))


  #Keep track of calendar months to be changed
  pmon <- which((par$mean != nochange_mean) | (par$var != nochange_var))

  # Non-zero days
  index_nz <- which(value>0)
  value_nz <- value[index_nz]
  mon.ts_nz <- mon.ts[index_nz]
  year.ts_nz <- year.ts[index_nz]

  # Non-zero days in each month
  value_nz_pmon <- lapply(1:12, function(x) value_nz[which(mon.ts_nz == x)])

  # Values to be changed (non-zero and ONLY perturbation months!)
  pind <- which((mon.ts %in% pmon) & (value > 0))

  #If no precip days, exit
  if (length(pind) == 0) return(value)

  # Fit base distribution and estimate parameters
  invisible(capture.output(bdist[["fit"]][pmon] <- lapply(pmon,
    function(x) fitdist(value_nz_pmon[[x]], "gamma"))))

  # Parameters for base distribution
  bdist[["shape"]][pmon] <- lapply(pmon, function(x) bdist[["fit"]][[x]]$estimate[[1]])
  bdist[["scale"]][pmon] <- lapply(pmon, function(x) 1/(bdist[["fit"]][[x]]$estimate[[2]]))
  bdist[["mean"]][pmon]  <- lapply(pmon, function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]])
  bdist[["var"]][pmon]   <- lapply(pmon, function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]]^2)

  # Define the mean, variance of the target distribution
  tdist[["mean"]][pmon, ] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) bdist[["mean"]][[m]] * parms_yearly[["mean"]][y, m]))
  tdist[["var"]][pmon, ]  <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) bdist[["var"]][[m]]  * parms_yearly[["var"]][y, m]))

  # Define the shape and scale parameters of the target distribution
  tdist[["scale"]][pmon,] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) tdist[["var"]][m,y] / tdist[["mean"]][m,y]))
  tdist[["shape"]][pmon,] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) tdist[["mean"]][m,y]^2 / tdist[["var"]][m,y]))

  # Estimate quantiles from base distribution
  qtile <- pgamma(value[pind], shape = unlist(bdist[["shape"]][mon.ts[pind]]), scale = unlist(bdist[["scale"]][mon.ts[pind]]))

  # Find time-series of shape and scale parameters for the fitted distribution
  shape_vec <- sapply(pind, function(x) tdist[["shape"]][mon.ts[x], year.ts[x]])
  scale_vec <- sapply(pind, function(x) tdist[["scale"]][mon.ts[x], year.ts[x]])

  # Replace the original matrix
  value[pind] <- qgamma(qtile, shape = shape_vec, scale = scale_vec)

  return(value)
}
