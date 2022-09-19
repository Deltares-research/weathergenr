
#' Function for perturbing climate statistics using Quantile-Mapping
#'
#' @param value To be completed...
#' @param mean.change placeholder
#' @param var.change placeholder
#' @param mon.ts placeholder
#' @param year.ts placeholder
#'
#' @return
#' @export
quantileMapping <- function(
  value = NULL,
  mean.change = NULL,
  var.change = NULL,
  mon.ts = NULL,
  year.ts = NULL,
  fit.method = "mme")

  {

  ymax <-  max(year.ts)
  emp1 <- matrix(0, nrow = 12, ncol = ymax)
  emp2 <- vector("list", 12)

  bdist <- list(fit = emp2, shape = emp2, scale = emp2, mean = emp2, var = emp2)
  tdist <- list(fit = NA, shape = emp1, scale = emp1, mean = emp1, var = emp1)

  #Keep track of calendar months to be changed
  pmon <- 1:12

  # Non-zero days
  value_nz <- value[which(value>0)]
  mon.ts_nz <- mon.ts[which(value>0)]

  # List of precip values in each month
  value_pmon <- lapply(1:12, function(x) value_nz[which(mon.ts_nz == x)])

  # Indices of events in months with at least 10 precip events
  pmon <- intersect(which(sapply(value_pmon, function(x) length(x) > 9)), pmon)

  # Values to be changed (non-zero and ONLY perturbation months!)
  pind <- which((mon.ts %in% pmon) & (value > 0))

  #If no precip days, exit
  if (length(pind) == 0) return(value)

  # Fit base distribution and estimate parameters
  bdist[["fit"]][pmon] <- lapply(pmon,
    function(x) fitdistrplus::fitdist(value_pmon[[x]], "gamma", method = fit.method))

  # Parameters for base distribution
  bdist[["shape"]][pmon] <- lapply(pmon,
    function(x) bdist[["fit"]][[x]]$estimate[[1]])
  bdist[["scale"]][pmon] <- lapply(pmon,
    function(x) 1/(bdist[["fit"]][[x]]$estimate[[2]]))
  bdist[["mean"]][pmon] <- lapply(pmon,
    function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]])
  bdist[["var"]][pmon] <- lapply(pmon,
    function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]]^2)

  # Parameters for the target distribution
  tdist[["mean"]][pmon, ] <- sapply(1:ymax, function(y)
    sapply(pmon, function(m) bdist[["mean"]][[m]] * mean.change[y, m]))
  tdist[["var"]][pmon, ] <- sapply(1:ymax, function(y)
    sapply(pmon, function(m) bdist[["var"]][[m]] * var.change[y, m]))

  # Define the shape and scale parameters of the target distribution
  tdist[["scale"]][pmon,] <- sapply(1:ymax, function(y)
    sapply(pmon, function(m) tdist[["var"]][m,y] / tdist[["mean"]][m,y]))
  tdist[["shape"]][pmon,] <- sapply(1:ymax, function(y)
    sapply(pmon, function(m) tdist[["mean"]][m,y]^2 / tdist[["var"]][m,y]))

  # Estimate quantiles from base distribution
  qtile <- stats::pgamma(value[pind],
    shape = unlist(bdist[["shape"]][mon.ts[pind]]),
    scale = unlist(bdist[["scale"]][mon.ts[pind]]))

  # Find time-series of shape and scale parameters for the fitted distribution
  shape_vec <- sapply(pind, function(x) tdist[["shape"]][mon.ts[x], year.ts[x]])
  scale_vec <- sapply(pind, function(x) tdist[["scale"]][mon.ts[x], year.ts[x]])

  # Replace the original matrix
  value[pind] <- stats::qgamma(qtile, shape = shape_vec, scale = scale_vec)

  return(value)
}
