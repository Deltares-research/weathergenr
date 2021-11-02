

#' A function for Wavelet Autoregressive Modeling (WARM)
#'
#' @param num.years A numeric value defining the desired length of simulated annual time-series
#' @param num.realizations A numeric value defining the number of synthetic series to be produced.
#' @param wavelet.components A list object, with different components corresponding to low-frequency signals and the noise
#'
#' @export
#' @import forecast
waveletARIMA <- function(
  wavelet.components = NULL,
  num.years = NULL,
  num.realizations = 1000)

  {

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

    SIM[[k]] <- replicate(num.realizations,
         simulate(MODEL[[k]], num.years, sd = sqrt(MODEL[[k]]$sigma2)) + INTERCEPT + MEAN)
  }

  return(Reduce(`+`, SIM))

}
