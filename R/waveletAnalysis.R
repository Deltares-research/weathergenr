
#' Wavelet Analysis Function
#'
#' \code{waveletAnalysis} returns the wavelet analysis results.
#'
#' To be completed later...
#'
#' @param variable              A numeric vector of time-series of variables, for example, time-series of annual precipitation.
#' @param signif.level          A numeric value to set the siginificance level of the wavelet analysis (Default= 0.90).
#' @param background.noise      A character string defining the type of background noise. Currently either "white" (default) or "red".
#' @param save.plot             A logical to save the results to file.
#' @param variable.unit         A character string to define the unit of the variable.
#'
#' @return
#' @export
#' @import ggplot2
#' @import cowplot
waveletAnalysis <- function(variable = NULL,
                            variable.unit = "mm",
                            signif.level = 0.90,
                            background.noise = "white",
                            plot = FALSE)
{


  #### Define Wavelet function used that returns transform
  waveletf <- function(k,s) {
    nn <- length(k)
    k0 <- 6    #nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
    z <- array(1,nn)
    z[which(k<=0)] <- 0
    expnt <- -((s*k - k0)^2/2)*z
    norm <- sqrt(s*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
    daughter <- norm*exp(expnt)
    daughter <- daughter*z
    fourier_factor <- (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
    coi <- fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
    dofmin <- 2  	# Degrees of freedom
    return(daughter)
  }

  #Define Wavelet function used that returns fourier_factor, cone of influence, and degrees of freedom
  waveletf2 <- function(k,s) {
    nn <- length(k)
    k0 <- 6    #nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
    z <- array(1,nn)
    z[which(k<=0)] <- 0
    expnt <- -((s*k - k0)^2/2)*z
    norm <- sqrt(s*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
    daughter <- norm*exp(expnt)
    daughter <- daughter*z
    fourier_factor <- (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
    coi <- fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
    dofmin <- 2  	# Degrees of freedom
    return(c(fourier_factor,coi,dofmin))
  }


  ######## PERFORM WAVELET  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #....construct time series to analyze, pad if necessary
  current_variable_org <- variable
  variance1 <- var(current_variable_org)
  n1 <- length(current_variable_org)

  current_variable <- scale(current_variable_org)
  variance2 <- var(current_variable)
  base2 <- floor(log(n1)/log(2) + 0.4999)   # power of 2 nearest to N

  current_variable <- c(current_variable,rep(0,(2^(base2+1)-n1)))
  n <- length(current_variable)

  #Determine parameters for Wavelet analysis
  dt <- 1
  dj <- 0.25
  s0 <- 2*dt
  J <- floor((1/dj)*log((n1*dt/s0),base=2))

  #....construct SCALE array & empty PERIOD & WAVE arrays
  scale <- s0*2^((0:J)*dj)
  period <- scale
  wave <- array(as.complex(0),c(J+1,n))  # define the wavelet array

  #....construct wavenumber array used in transform [Eqn(5)]
  k <- c(1:floor(n/2))
  k <- k*((2.*pi)/(n*dt))
  k <- c(0,k,-rev(k[1:floor((n-1)/2)]))

  #fourier transform of standardized precipitation
  f <- fft(current_variable,inverse=FALSE)

  # loop through all scales and compute transform
  for (a1 in 1:(J+1)) {
    daughter <- waveletf(k,scale[a1])
    results <- waveletf2(k,scale[a1])
    fourier_factor <- results[1]
    coi <- results[2]
    dofmin <- results[3]
    wave[a1,] <- fft(f*daughter,inverse=TRUE)/n  # wavelet transform[Eqn(4)]
  }

  period <- fourier_factor*scale
  coi <- coi*dt*c((0.00001),1:((n1+1)/2-1),rev((1:(n1/2-1))),(0.00001))  # COI [Sec.3g]
  wave <- wave[,1:n1]  # get rid of padding before returning
  POWER <- abs(wave)^2
  GWS <- variance1*apply(POWER,FUN=mean,c(1)) #Global Wavelet Spectrum


  ######## Signficance Testing  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  # get the appropriate parameters [see Table(2)]
  k0 <- 6
  empir <- c(2,0.776,2.32,0.60)
  dofmin <- empir[1]     # Degrees of freedom with no smoothing
  Cdelta <- empir[2]     # reconstruction factor
  gamma_fac <- empir[3]  # time-decorrelation factor
  dj0 <- empir[4]       # scale-decorrelation factor

  if (background.noise=="white") {lag1 <- 0}
  if (background.noise=="red") {lag1 <- .72}

  freq <- dt / period   # normalized frequency
  fft_theor <- (1-lag1^2) / (1-2*lag1*cos(freq*2*pi)+lag1^2)  # [Eqn(16)]
  fft_theor <- fft_theor  # include time-series variance
  dof <- dofmin

  #ENTIRE POWER SPECTRUM
  chisquare <- qchisq(signif.level,dof)/dof
  signif <- fft_theor*chisquare   # [Eqn(18)]
  sig95 <- ((signif))%o%(array(1,n1))  # expand signif --> (J+1)x(N) array
  sig95 <- POWER / sig95         # where ratio > 1, power is significant

  #TIME_AVERAGED (GLOBAL WAVELET SPECTRUM)
  dof <- n1 - scale
  if (length(dof) == 1) {dof <- array(0,(J+1))+dof}
  dof[which(dof < 1)] <- 1
  dof <- dofmin*sqrt(1 + (dof*dt/gamma_fac / scale)^2 )
  tt <- which(dof < dofmin)
  dof[tt] <- dofmin
  chisquare_GWS <- array(NA,(J+1))
  GWS_signif <- array(NA,(J+1))

  for (a1 in 1:(J+1)) {
    chisquare_GWS[a1] <- qchisq(signif.level,dof[a1])/dof[a1]
    GWS_signif[a1] <- fft_theor[a1]*variance1*chisquare_GWS[a1]
  }

  period_lower_limit <- 0
  sig_periods <- which(GWS>GWS_signif & period > period_lower_limit)
  signif_periods <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))

  if(plot == TRUE) {

    GWS_gg <- tibble(period = period, GWS = GWS, GWS_signif = GWS_signif)
    var_gg <- tibble(x = 1:length(variable), y = variable)


    df <- t(log(POWER,base=2)) %>%
      as_tibble(.name_repair = ~ as.character(period)) %>%
      mutate(x = 1:n1) %>%
      gather(key = y, value = z, -x) %>%
      mutate(across(everything(), as.numeric))

    df <- with(df,
               akima::interp(x, y, z, extrap = FALSE, linear = TRUE,
                      xo = seq(min(x), max(x), length = 20),
                      yo = seq(min(y), max(y), length = 20)))

    df1 <- as_tibble(df$z)  %>%
      setNames(df$y) %>% mutate(x = df$x) %>%
      gather(key = y, value = z, -x) %>%
      mutate(across(everything(), as.numeric))


    df2 <- tibble(x= 1:n1, y = coi) %>% filter(y > min(df1$x))

    df3 <- as_tibble(t(sig95)) %>%
      setNames(period) %>% mutate(x = 1:n1) %>%
      gather(key = y, value = z, -x) %>%
      mutate(across(everything(), as.numeric))


    p1 <- ggplot(var_gg) +
      theme_light() +
      geom_line(aes(x,y)) +
      geom_point(aes(x,y), shape = 21, size = 2) +
      labs(x= "Time (years)", y = variable.unit) +
      scale_x_continuous() +
      ggtitle("Annual Time-series")

    p2 <- ggplot(df1, aes(x=x,y=y)) +
      theme_light() +
      geom_raster(aes(fill = z)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_reverse(expand=c(0,0)) +
      viridis::scale_fill_viridis() +
      labs(x= "Time (years)", y = "Period (years)") +
      guides(fill = "none") +
      geom_line(data = df2, linetype = "dashed") +
      stat_contour(aes(z = z), data = df3, breaks=c(-99, 1)) +
      ggtitle("Wavelet Power Spectrum")

    p3 <- ggplot(GWS_gg) +
      theme_light() +
      geom_line(aes(period, GWS)) +
      geom_point(aes(period, GWS), shape = 21, size = 2) +
      geom_line(aes(period, GWS_signif), color = "red", linetype = "dashed") +
      scale_y_continuous() +
      scale_x_reverse(expand = c(0,0)) +
      coord_flip() +
      labs(y = expression(Power~(mm^2)), x="") +
      ggtitle("Global Wavelet Power Spectrum")


    p23 <- plot_grid(p2, p3, nrow = 1, rel_widths = c(2,1), labels = c("b)", "c)"))
    p <- plot_grid(p1, p23, nrow = 2, rel_heights = c(1,2.5), labels = c("a)", ""))

  } else {
    p <- NULL
  }

  return(list(GWS = GWS,
              GWS_signif = GWS_signif,
              GWS_period = period,
              signif_periods = signif_periods,
              plot = p))

}

