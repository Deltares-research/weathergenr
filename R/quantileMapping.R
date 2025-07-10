#' Quantile Mapping for Climate Change Perturbation of Precipitation
#'
#' @description
#' Applies quantile mapping to adjust a daily precipitation time series to reflect changes in monthly mean and variance,
#' consistent with climate change scenarios. The method fits a gamma distribution to observed nonzero precipitation
#' in each calendar month, then modifies the distribution's mean and variance according to supplied monthly/yearly change factors,
#' and maps the original values to their new quantiles in the perturbed distribution.
#'
#' This function is designed to be used within larger climate stress-testing workflows to impose user-specified
#' changes in mean and variance on daily precipitation data, while preserving realistic distributional characteristics.
#'
#' @param value Numeric vector. Original daily precipitation values to perturb (typically for one grid cell), must be non-negative.
#' @param mean.change Matrix or array of mean change factors (multiplicative), dimension: `n_years` x 12 (year, month). Each entry indicates the scaling for that year's/month's mean.
#' @param var.change Matrix or array of variance change factors (multiplicative), dimension: `n_years` x 12 (year, month).
#' @param mon.ts Integer vector (same length as `value`). Calendar month for each day (1-12).
#' @param year.ts Integer vector (same length as `value`). Simulation year index for each day (1 = first year, etc).
#' @param fit.method Character. Method for fitting the base gamma distribution; passed to [fitdistrplus::fitdist()]. Default is `"mme"`.
#'
#' @return
#' Numeric vector, same length as `value`. Precipitation time series perturbed according to quantile mapping procedure.
#'
#' @details
#' - Only nonzero precipitation days in months with at least 10 such days are perturbed; all others are returned unchanged.
#' - For each month and simulation year, a gamma distribution is fitted to the base data, then mean and variance are scaled
#'   by `mean.change` and `var.change` to define the perturbed distribution.
#' - The original nonzero value is mapped to a quantile in the base distribution, then transformed to the same quantile
#'   in the target distribution.
#' - The `mean.change` and `var.change` inputs should be either matrices with dimensions (`n_years` x 12), or objects
#'   that can be indexed as `mean.change[year, month]`.
#'
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pgamma qgamma
#'
#' @examples
#' \dontrun{
#' # Example: 2 years of daily data, simple 5% mean and 10% variance increase
#' set.seed(123)
#' n_days <- 730
#' year_idx <- rep(1:2, each = 365)
#' month_idx <- rep(rep(1:12, times = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)), 2)[1:n_days]
#' daily_precip <- rgamma(n_days, shape = 1, scale = 5)
#'
#' mean.change <- matrix(1.05, nrow = 2, ncol = 12)
#' var.change <- matrix(1.10, nrow = 2, ncol = 12)
#'
#' perturbed_precip <- quantileMapping(
#'   value = daily_precip,
#'   mean.change = mean.change,
#'   var.change = var.change,
#'   mon.ts = month_idx,
#'   year.ts = year_idx
#' )
#'
#' summary(daily_precip)
#' summary(perturbed_precip)
#' }
#'
#' @seealso \code{\link[fitdistrplus]{fitdist}}
#' @export
quantileMapping <- function(
    value = NULL,
    mean.change = NULL,
    var.change = NULL,
    mon.ts = NULL,
    year.ts = NULL,
    fit.method = "mme") {
  ymax <- max(year.ts)
  emp1 <- matrix(0, nrow = 12, ncol = ymax)
  emp2 <- vector("list", 12)

  bdist <- list(fit = emp2, shape = emp2, scale = emp2, mean = emp2, var = emp2)
  tdist <- list(fit = NA, shape = emp1, scale = emp1, mean = emp1, var = emp1)

  # Keep track of calendar months to be changed
  pmon <- 1:12

  # Non-zero days
  value_nz <- value[which(value > 0)]
  mon.ts_nz <- mon.ts[which(value > 0)]

  # List of precip values in each month
  value_pmon <- lapply(1:12, function(x) value_nz[which(mon.ts_nz == x)])

  # Indices of events in months with at least 10 precip events
  pmon <- intersect(which(sapply(value_pmon, function(x) length(x) > 9)), pmon)

  # Values to be changed (non-zero and ONLY perturbation months!)
  pind <- which((mon.ts %in% pmon) & (value > 0))

  # If no precip days, exit
  if (length(pind) == 0) {
    return(value)
  }

  # Fit base distribution and estimate parameters
  bdist[["fit"]][pmon] <- lapply(
    pmon,
    function(x) fitdistrplus::fitdist(value_pmon[[x]], "gamma", method = fit.method)
  )

  # Parameters for base distribution
  bdist[["shape"]][pmon] <- lapply(
    pmon,
    function(x) bdist[["fit"]][[x]]$estimate[[1]]
  )
  bdist[["scale"]][pmon] <- lapply(
    pmon,
    function(x) 1 / (bdist[["fit"]][[x]]$estimate[[2]])
  )
  bdist[["mean"]][pmon] <- lapply(
    pmon,
    function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]]
  )
  bdist[["var"]][pmon] <- lapply(
    pmon,
    function(x) bdist[["shape"]][[x]] * bdist[["scale"]][[x]]^2
  )

  # Parameters for the target distribution
  tdist[["mean"]][pmon, ] <- sapply(1:ymax, function(y) {
    sapply(pmon, function(m) bdist[["mean"]][[m]] * mean.change[y, m])
  })
  tdist[["var"]][pmon, ] <- sapply(1:ymax, function(y) {
    sapply(pmon, function(m) bdist[["var"]][[m]] * var.change[y, m])
  })

  # Define the shape and scale parameters of the target distribution
  tdist[["scale"]][pmon, ] <- sapply(1:ymax, function(y) {
    sapply(pmon, function(m) tdist[["var"]][m, y] / tdist[["mean"]][m, y])
  })
  tdist[["shape"]][pmon, ] <- sapply(1:ymax, function(y) {
    sapply(pmon, function(m) tdist[["mean"]][m, y]^2 / tdist[["var"]][m, y])
  })

  # Estimate quantiles from base distribution
  qtile <- stats::pgamma(value[pind],
    shape = unlist(bdist[["shape"]][mon.ts[pind]]),
    scale = unlist(bdist[["scale"]][mon.ts[pind]])
  )

  # Find time-series of shape and scale parameters for the fitted distribution
  shape_vec <- sapply(pind, function(x) tdist[["shape"]][mon.ts[x], year.ts[x]])
  scale_vec <- sapply(pind, function(x) tdist[["scale"]][mon.ts[x], year.ts[x]])

  # Replace the original matrix
  value[pind] <- stats::qgamma(qtile, shape = shape_vec, scale = scale_vec)

  return(value)
}
