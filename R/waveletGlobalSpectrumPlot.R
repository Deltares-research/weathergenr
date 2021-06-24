

#' @title ggplot wrapper for wavelet spectra plot
#'
#' @description This function is simply a ggplot wrapper for displaying results from a wavelet analysis.
#'
#' @param power.period A numeric vector of fourier periods.
#' @param power.signif A numeric value to set the sigicance level.
#' @param power.obs A numeric vector of observed global wavelet spectra.
#' @param power.sim A numeric matrix of simulated global wavelet spectrum.
#'
#' @return A ggplot2 object
#' @details place holder
#' @export

waveletGlobalSpectrumPlot <- function(power.period, power.signif, power.obs, power.sim = NULL) {

  require(ggplot2)
  require(scales)
  require(dplyr)

  tsn <- length(power.period)

  if(is.null(power.sim)) {

    df <- tibble(power.period, power.signif, power.obs)

    p <- ggplot(df, aes(x = power.period)) +
      geom_line(aes(y=power.signif), color = "red", linetype = "dashed", size=0.6) +
      geom_line(aes(y=power.obs), color = "blue", size=0.6)

  } else {

    savg <- apply(power.sim, 1, mean)
    slow <- apply(power.sim, 1, min)
    sup  <- apply(power.sim, 1, max)

    #slow <- apply(power.sim, 1, function(x) quantile(x, 0.025))
    #sup  <- apply(power.sim, 1, function(x) quantile(x, 0.975))

    df <- tibble(power.period, power.signif, power.obs, slow=slow[1:tsn], sup=sup[1:tsn], savg=savg[1:tsn])


    p <- ggplot(df, aes( x = power.period)) +
      geom_ribbon(aes(ymin = slow, ymax = sup), alpha = 0.2) +
      geom_line(aes(y=power.signif), color = "red" , linetype = "dashed", size=0.6) +
      geom_line(aes(y=savg), color = "black", size=0.6) +
      geom_line(aes(y=power.obs), color = "blue", size=0.6)
  }

  p <- p + theme_light(base_size = 11) +
    labs(x="Period (years)", y = expression(paste("Power (", mm^2,")"))) #+
    #scale_x_continuous(breaks=seq(0,50,5), expand=c(0,0)) +
    #scale_y_log10(limits = c(5*10**2,5*10**5),labels = comma)

  return(p)

}
