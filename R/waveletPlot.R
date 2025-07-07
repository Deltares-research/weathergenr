

#' @title A ggplot wrapper for wavelet spectra plot
#'
#' @description ggplot wrapper for visualization of wavelet spectral analysis
#'
#' @param power.period A numeric vector of fourier periods from wavelet analysis
#' @param power.signif A numeric value to set the significance level of the spectral analysis
#' @param power.obs A numeric vector of observed  global wavelet spectra.
#' @param power.sim A numeric matrix of simulated global wavelet spectrum. Each column is an independent observation
#'
#' @return A ggplot2 object
#' @details place holder
#' @export
#' @import ggplot2 scales dplyr
waveletPlot <- function(
  power.period = NULL,
  power.signif = NULL,
  power.obs = NULL,
  power.sim = NULL)
  {

  tsn <- length(power.period)

  if(is.null(power.sim)) {

    df <- tibble(power.period, power.signif, power.obs)

    p <- ggplot(df, aes(x = power.period)) +
      theme_light(base_size = 11) +
      geom_line(aes(y=power.signif), color = "red", linetype = "dashed", linewidth = 0.6) +
      geom_line(aes(y=power.obs), color = "blue", linewidth = 0.6)

  } else {

    savg <- apply(as.matrix(power.sim), 1, mean)
    slow <- apply(as.matrix(power.sim), 1, min)
    sup  <- apply(as.matrix(power.sim), 1, max)

    #slow <- apply(power.sim, 1, function(x) quantile(x, 0.025))
    #sup  <- apply(power.sim, 1, function(x) quantile(x, 0.975))

    df <- tibble(power.period, power.signif, power.obs, slow=slow[1:tsn], sup=sup[1:tsn], savg=savg[1:tsn])

    p <- ggplot(df, aes( x = power.period)) +
      theme_light(base_size = 11) +
      geom_ribbon(aes(ymin = slow, ymax = sup), alpha = 0.2) +
      geom_line(aes(y=power.signif), color = "red" , linetype = "dashed", linewidth = 0.6) +
      geom_line(aes(y=savg), color = "black", linewidth = 0.6) +
      geom_line(aes(y=power.obs), color = "blue", linewidth = 0.6)
  }

  p <- p + theme_light(base_size = 11) +
    labs(x="Period (years)", y = expression(paste("Power (", mm^2,")"))) #+
    #scale_x_continuous(breaks=seq(0,50,5), expand=c(0,0)) +
    #scale_y_log10(limits = c(5*10**2,5*10**5),labels = comma)

  return(p)

}
