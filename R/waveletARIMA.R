#' Wavelet Autoregressive Modeling (WARM)
#'
#' Simulates synthetic time series by modeling each wavelet component (signal or noise)
#' with an ARIMA model, and then summing the simulated components.
#'
#' @param wavelet.components A list or matrix where each column (or list element) is a numeric vector corresponding to a wavelet component (low-frequency signal or noise).
#' @param sim.year.num Integer. Desired length (number of years or timesteps) of each simulated series.
#' @param sim.num Integer. Number of synthetic series to produce. Default: 1000.
#' @param seed Optional. Integer random seed for reproducibility.
#'
#' @return A matrix of dimension \code{sim.year.num} x \code{sim.num}, where each column is a synthetic time series realization.
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate 10 synthetic annual series (length 30) from two "wavelet" components
#' set.seed(42)
#' wavelet.components <- list(
#'   signal = sin(2 * pi * 1:30 / 8) + rnorm(30, 0, 0.2), # Low-freq
#'   noise  = rnorm(30, 0, 0.5)
#' )
#' result <- waveletARIMA(wavelet.components, sim.year.num = 30, sim.num = 10, seed = 123)
#' dim(result)  # Should be 30 x 10
#' matplot(result, type = "l", lty = 1, main = "Synthetic Wavelet-ARIMA Realizations")
#' }
waveletARIMA <- function(
    wavelet.components = NULL,
    sim.year.num = NULL,
    sim.num = 1000,
    seed = NULL
) {

  # Input validation
  if (is.null(wavelet.components)) {
    stop("Input 'wavelet.components' must not be NULL.")
  }
  if (is.data.frame(wavelet.components) || is.matrix(wavelet.components)) {
    # For data.frames, check for numeric columns only
    if (any(!vapply(wavelet.components, is.numeric, logical(1)))) {
      stop("All columns of 'wavelet.components' must be numeric.")
      return(invisible(NULL)) # Ensure function stops
    }
  } else if (is.list(wavelet.components)) {
    if (any(!vapply(wavelet.components, is.numeric, logical(1)))) {
      stop("All elements of 'wavelet.components' must be numeric vectors.")
      return(invisible(NULL))
    }
  } else {
    stop("Input 'wavelet.components' must be a matrix, data.frame, or list of numeric vectors.")
    return(invisible(NULL))
  }

  # Accept both list or matrix/data.frame
  if (is.matrix(wavelet.components) || is.data.frame(wavelet.components)) {
    ncomp <- ncol(wavelet.components)
    comp_list <- lapply(seq_len(ncomp), function(k) wavelet.components[, k])
  } else if (is.list(wavelet.components)) {
    ncomp <- length(wavelet.components)
    comp_list <- wavelet.components
  } else {
    stop("'wavelet.components' must be a list, matrix, or data.frame")
  }

  # Simulate each component
  SIM <- vector("list", ncomp)
  for (k in seq_len(ncomp)) {

    component <- unlist(comp_list[[k]], use.names = FALSE)
    MEAN <- mean(component)
    CENTERED <- component - MEAN

    # Fit ARIMA with auto.arima, restrict to non-seasonal ARMA
    suppressPackageStartupMessages({
      MODEL <- forecast::auto.arima(CENTERED,
        max.p = 2, max.q = 2, max.P = 0, max.Q = 0,
        stationary = TRUE, seasonal = FALSE)
    })
    INTERCEPT <- if ("intercept" %in% names(MODEL$coef)) MODEL$coef[["intercept"]] else 0
    SD <- sqrt(MODEL$sigma2)

    # Simulate sim.num synthetic series for this component
    if (!is.null(seed)) {
      old_seed <- .Random.seed
      on.exit({ .Random.seed <<- old_seed }, add = TRUE)
      set.seed(seed + k*1000)
    }

    SIM[[k]] <- replicate(sim.num, {
      # To ensure different random streams, combine seed and simulation index
      as.numeric(stats::simulate(MODEL, sim.year.num, sd = SD)) + INTERCEPT + MEAN
    })

}

  # Sum across components: each SIM[[k]] is sim.year.num x sim.num
  if (ncomp == 1) {
    output <- SIM[[1]]
  } else {
    output <- Reduce(`+`, SIM)
  }
  # Return as a matrix: sim.year.num x sim.num
  return(output)
}

