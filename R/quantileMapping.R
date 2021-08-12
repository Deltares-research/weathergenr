
#' Function for perturbing climate statistics using Quantile-Mapping
#'
#' @param par List of parameters to be passed.
#' @param operator A character string to specify the type of perturb change operator.
#' @param value To be completed...
#' @param date  To be completed...
#'
#' @return
#' @export
#'
quantileMapping <- function(
  value = NULL,
  date = NULL,
  par = NULL,
  operator = "multiply")

  {

  #Date objects
  mon.ts <- month(date)
  year.ts <- year(date) - min(year(date)) + 1

  emp <- matrix(0, nrow = 12, ncol = max(year.ts))

  bdist <- list(fit = NA, shape = NA, scale = NA, mean = NA, var = NA)
  tdist <- list(fit = NA, shape = emp, scale = emp, mean = emp, var = emp)

  # Identify months to be perturbed
  nochange_var = 1
  if(operator == "multiply") {
    nochange_mean <- 1
  } else if(operator == "add") {
    nochange_mean <- 0
  } else {
    stop("change operator can be either 'multiply' or 'add'")
  }

  #Keep track of changing months
  pmon <- which((par$mean != nochange_mean) | (par$var != nochange_var))

  # Values to be changes (non-zero days for only perturbed months!)
  pind <- which((mon.ts %in% pmon) & (value > 0))

  #If no precip days, exit
  if (length(pind) == 0) {
    return(value)

  } else {

    pind_permon <- sapply(1:12, function(x) which(month(date[pind]) == x))

    # Fit base distribution and estimate parameters
    bdist[["fit"]][pmon] <- lapply(pmon, function(x) quiet(fitdist(value[pind][pind_permon[[x]]], distr = "gamma", method = "mle")))
    bdist[["shape"]][pmon] <- sapply(pmon, function(x) bdist[["fit"]][[x]]$estimate[[1]])
    bdist[["scale"]][pmon] <- sapply(pmon, function(x) 1/(bdist[["fit"]][[x]]$estimate[[2]]))
    bdist[["mean"]][pmon] <- sapply(pmon, function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]])
    bdist[["var"]][pmon] <- sapply(pmon, function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]]^2)

    # Vector of means and variances (row=years, column=months)
    parms_yearly <- list()
    parms_yearly[["mean"]] <- sapply(1:12, function(m) seq(nochange_mean, par$mean[m], length.out = max(year.ts)))
    parms_yearly[["var"]] <- sapply(1:12, function(m) seq(nochange_var, par$var[m], length.out = max(year.ts)))

    if (operator == "multiply") {

      # Define the mean, variance of the target distribution
      tdist[["mean"]][pmon, ] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) bdist[["mean"]][m] * parms_yearly[["mean"]][y, m]))
      tdist[["var"]][pmon, ] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) bdist[["var"]][m] * parms_yearly[["var"]][y, m]))

    } else {

      # Define the mean, variance of the target distribution
      tdist[["mean"]][pmon,] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) bdist[["mean"]][[m]] + parms_yearly[["mean"]][y, m]))
      tdist[["var"]][pmon,]  <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) bdist[["var"]][[m]] * parms_yearly[["var"]][y, m]))

    }

    # Define the shape and scale parameters of the target distribution
    tdist[["scale"]][pmon,] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) tdist[["var"]][m,y] / tdist[["mean"]][m,y]))
    tdist[["shape"]][pmon,] <- sapply(1:max(year.ts), function(y) sapply(pmon, function(m) tdist[["mean"]][m,y]^2 / tdist[["var"]][m,y]))


    # Estimate quantiles from base distribution
    qtile <- pgamma(value[pind], shape = bdist[["shape"]][mon.ts[pind]], scale = bdist[["scale"]][mon.ts[pind]])

    # Find time-series of shape and scale parameters for the fitted distribution
    shape_vec <- sapply(pind, function(x) tdist[["shape"]][mon.ts[x], year.ts[x]])
    scale_vec <- sapply(pind, function(x) tdist[["scale"]][mon.ts[x], year.ts[x]])

    # Replace the original matrix
    value[pind] <- qgamma(qtile, shape = shape_vec, scale = scale_vec)

    return(value)

  }
}

  # # Empty lists and vectors to store distributions
  # emp1 <- vector("list", length = 12)
  # bdist <- list(fit = emp1, shape = emp1, scale = emp1, mean = emp1, var = emp1)
  #
  # emp2 <- matrix(0, nrow = 12, ncol = max(year.ts))
  # tdist <- list(fit = emp2, shape = emp2, scale = emp2, mean = emp2, var = emp2)
