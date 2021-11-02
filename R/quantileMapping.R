
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
quantileMapping <- function(
  value = NULL,
  date = NULL,
  par = NULL,
  mon.ts = NULL,
  year.ts = NULL,
  step.change = TRUE)

  {

  ymax <-  max(year.ts)
  emp1 <- matrix(0, nrow = 12, ncol = ymax)
  emp2 <- vector("list", 12)

  bdist <- list(fit = emp2, shape = emp2, scale = emp2, mean = emp2, var = emp2)
  tdist <- list(fit = NA, shape = emp1, scale = emp1, mean = emp1, var = emp1)

  # Vector of means and variances (row=years, column=months)
  par_a <- list()

  if(step.change) {
    par_a[["mean"]] <- sapply(1:12, function(m) seq(1, par$mean[m], length.out = ymax))
    par_a[["var"]]  <- sapply(1:12, function(m) seq(1, par$var[m], length.out = ymax))
  } else {
    par_a[["mean"]] <- sapply(1:12, function(m) rep(par$mean[m], ymax))
    par_a[["var"]]  <- sapply(1:12, function(m) rep(par$var[m], ymax))
  }

  #Keep track of calendar months to be changed
  pmon <- which((par$mean != 1) | (par$var != 1))

  # Non-zero days
  index_nz <- which(value>0)
  value_nz <- value[index_nz]
  mon.ts_nz <- mon.ts[index_nz]
  year.ts_nz <- year.ts[index_nz]

  # Non-zero days in each month
  value_pmon <- lapply(1:12, function(x) value_nz[which(mon.ts_nz == x)])

  # Values to be changed (non-zero and ONLY perturbation months!)
  pind <- which((mon.ts %in% pmon) & (value > 0))

  #If no precip days, exit
  if (length(pind) == 0) return(value)

  # Fit base distribution and estimate parameters
  bdist[["fit"]][pmon] <- lapply(pmon,
    function(x) fitdistrplus::fitdist(value_pmon[[x]], "gamma", control = list(reltol=1e-3)))

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
    sapply(pmon, function(m) bdist[["mean"]][[m]] * par_a[["mean"]][y, m]))
  tdist[["var"]][pmon, ] <- sapply(1:ymax, function(y)
    sapply(pmon, function(m) bdist[["var"]][[m]] * par_a[["var"]][y, m]))

  # Define the shape and scale parameters of the target distribution
  tdist[["scale"]][pmon,] <- sapply(1:ymax, function(y)
    sapply(pmon, function(m) tdist[["var"]][m,y] / tdist[["mean"]][m,y]))
  tdist[["shape"]][pmon,] <- sapply(1:ymax, function(y)
    sapply(pmon, function(m) tdist[["mean"]][m,y]^2 / tdist[["var"]][m,y]))

  # Estimate quantiles from base distribution
  qtile <- stats::pgamma(value[pind], shape = unlist(bdist[["shape"]][mon.ts[pind]]),
    scale = unlist(bdist[["scale"]][mon.ts[pind]]))

  # Find time-series of shape and scale parameters for the fitted distribution
  shape_vec <- sapply(pind, function(x) tdist[["shape"]][mon.ts[x], year.ts[x]])
  scale_vec <- sapply(pind, function(x) tdist[["scale"]][mon.ts[x], year.ts[x]])

  # Replace the original matrix
  value[pind] <- stats::qgamma(qtile, shape = shape_vec, scale = scale_vec)

  return(value)
}
