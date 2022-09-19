

#' Wavelet Autoregressive Modeling (WARM)
#'
#' @param sim.year.num A numeric value defining the desired length of simulated annual time-series
#' @param sim.num A numeric value defining the number of synthetic series to be produced.
#' @param wavelet.components A list object, with different components corresponding to low-frequency signals and the noise
#' @param seed random seed value
#'
#' @export
waveletARIMA <- function(
  wavelet.components = NULL,
  sim.year.num = NULL,
  sim.num = 1000,
  seed = NULL)
  {

  # If ne seed provided, generate a random number
  if(is.null(seed)) seed = sample.int(1e10,1)

  # Define ARIMA model for each component
  MODEL <- vector(mode = "list", length = ncol(wavelet.components))
  SIM   <- vector(mode = "list", length = ncol(wavelet.components))

  for (k in 1:ncol(wavelet.components)) {

    # Remove the mean from the component
    MEAN  <- mean(wavelet.components[[k]])
    CENTERED  <- wavelet.components[[k]] - MEAN

    MODEL[[k]] <-  forecast::auto.arima(CENTERED, max.p = 2,max.q = 2,max.P = 0,max.Q = 0,
      stationary = TRUE, seasonal = FALSE)

    INTERCEPT <- ifelse(length(which(names(MODEL[[k]]$coef)=="intercept")) > 0,
                        as.vector(MODEL[[k]]$coef)[which(names(MODEL[[k]]$coef)=="intercept")],0)

    SD <- sqrt(MODEL[[k]]$sigma2)

    SIM[[k]] <- sapply(1:sim.num, function(x) {
      set.seed(seed+x)
      stats::simulate(MODEL[[k]], sim.year.num, sd = SD)}) + INTERCEPT + MEAN
  }

  return(Reduce(`+`, SIM))

  }
