#' Wavelet Spectral Analysis with Significance Testing
#'
#' Performs continuous wavelet transform analysis on a time series,
#' computes the global wavelet spectrum (GWS) and significance levels,
#' and (optionally) produces a plot of power spectrum.
#'
#' @param variable Numeric vector, the time series (e.g. annual precipitation).
#' @param signif.level Numeric between 0 and 1. Significance level for test (default = 0.90).
#' @param noise.type "white" (default) or "red". Background noise for spectrum test.
#' @param variable.unit Character; unit for variable (for plotting).
#' @param plot Logical, whether to save a plot (default FALSE).
#' @param output.path Optional string. Directory where plot will be saved (if plot = TRUE).
#'
#' @return A list with:
#'   \item{GWS}{Global Wavelet Spectrum}
#'   \item{GWS_signif}{Significance threshold for GWS}
#'   \item{GWS_period}{Period axis for GWS}
#'   \item{signif_periods}{Indices of significant periods}
#'
#' @examples
#' set.seed(123)
#' # Simulate AR(1) process with periodicity at 8 years
#' years <- 1:64
#' signal <- sin(2 * pi * years / 8) + arima.sim(n=64, model=list(ar=0.7))
#' res <- waveletAnalysis(signal, plot=FALSE)
#' print(res$GWS_period[res$signif_periods])
#'
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @export
waveletAnalysis <- function(variable,
                            signif.level = 0.90,
                            noise.type = "white",
                            variable.unit = "mm",
                            plot = FALSE,
                            output.path = NULL)
{


  # --- Input checks
  stopifnot(is.numeric(variable), length(variable) > 8)
  if (anyNA(variable)) stop("Variable contains missing values.")
  if (!(noise.type %in% c("white", "red"))) stop("noise.type must be 'white' or 'red'")
  if (!is.numeric(signif.level) || signif.level <= 0 || signif.level >= 1)
    stop("signif.level must be between 0 and 1.")

  # Workaround for rlang warning
  x <- y <- z <- 0

  # Morlet wavelet
  waveletf <- function(k, s) {
    nn <- length(k)
    k0 <- 6
    z <- as.numeric(k > 0)
    expnt <- -((s * k - k0)^2 / 2) * z
    norm <- sqrt(s * k[2]) * (pi^(-0.25)) * sqrt(nn)
    daughter <- norm * exp(expnt) * z
    return(daughter)
  }

  waveletf2 <- function(k, s) {
    k0 <- 6
    fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0^2))
    coi <- fourier_factor / sqrt(2)
    dofmin <- 2
    c(fourier_factor, coi, dofmin)
  }

  # Standardize & pad
  variable_org <- variable
  variance1 <- stats::var(variable_org)
  n1 <- length(variable_org)
  variable <- scale(variable_org)
  base2 <- floor(log2(n1) + 0.4999)
  variable <- c(variable, rep(0, (2^(base2 + 1) - n1)))
  n <- length(variable)

  # Wavelet transform parameters
  dt <- 1; dj <- 0.25; s0 <- 2 * dt
  J <- floor((1/dj) * log((n1 * dt / s0), base=2))
  scale <- s0 * 2^((0:J) * dj)
  k <- c(0:(floor(n/2)), -rev(1:floor((n-1)/2))) * ((2 * pi) / (n * dt))


  # FFT of series
  f <- stats::fft(variable)
  wave <- array(as.complex(0), c(J+1, n))
  for (a1 in 1:(J+1)) {
    daughter <- waveletf(k, scale[a1])
    wave[a1, ] <- stats::fft(f * daughter, inverse = TRUE) / n
    if (a1 == 1) {
      params <- waveletf2(k, scale[a1])
      fourier_factor <- params[1]; coi_base <- params[2]
    }
  }
  period <- fourier_factor * scale
  coi <- coi_base * dt * c(0.00001, 1:((n1+1)/2-1), rev(1:(n1/2-1)), 0.00001)
  wave <- wave[, 1:n1]
  POWER <- abs(wave)^2
  GWS <- variance1 * rowMeans(POWER)

  # --- Significance testing
  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1]; gamma_fac <- empir[3]
  lag1 <- if (noise.type == "white") 0 else 0.72
  freq <- dt / period
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2)
  chisquare <- stats::qchisq(signif.level, dofmin) / dofmin
  signif <- fft_theor * chisquare
  sig95 <- POWER / (outer(signif, rep(1, n1)))
  dof <- n1 - scale
  dof[dof < 1] <- 1
  dof <- dofmin * sqrt(1 + (dof * dt / gamma_fac / scale)^2)
  dof[dof < dofmin] <- dofmin
  chisquare_GWS <- stats::qchisq(signif.level, dof) / dof
  GWS_signif <- fft_theor * variance1 * chisquare_GWS

  period_lower_limit <- 0
  sig_periods <- which(GWS > GWS_signif & period > period_lower_limit)
  sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))
  signif_periods <- unlist(lapply(sig_periods_grp, function(x) as.integer(round(median(x)))), use.names = FALSE)

  # --- Plotting (optional)
  if (plot) {

    if (is.null(output.path)) output.path <- tempdir()
    if (!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)

    GWS_gg <- data.frame(period = period, GWS = GWS, GWS_signif = GWS_signif)
    var_gg <- data.frame(x = 1:length(variable), y = variable)

    df <- t(log(POWER,base=2)) %>%
      tibble::as_tibble(.name_repair = ~ as.character(period)) %>%
      dplyr::mutate(x = 1:n1) %>%
      tidyr::gather(key = y, value = z, -x) %>%
      dplyr::mutate(across(everything(), as.numeric))

    df <- suppressWarnings(with(df, akima::interp(x, y, z, extrap = TRUE, linear = FALSE,
                                                  xo = seq(min(x), max(x), length = 20),
                                                  yo = seq(min(y), max(y), length = 20))))


    df1 <- tibble::as_tibble(df$z, .name_repair = ~as.character(df$y))  %>%
      dplyr::mutate(x = df$x) %>%
      tidyr::gather(key = y, value = z, -x) %>%
      dplyr::mutate(across(everything(), as.numeric))


    df2 <- tibble::tibble(x= 1:n1, y = coi) %>% filter(y > min(df1$x))

    df3 <- tibble::as_tibble(t(sig95), .name_repair = ~as.character(period)) %>%
      dplyr::mutate(x = 1:n1) %>%
      tidyr::gather(key = y, value = z, -x) %>%
      dplyr::mutate(across(everything(), as.numeric))

    p2 <- ggplot2::ggplot(df1, ggplot2::aes(x = x, y = y)) +
      ggplot2::theme_light() +
      ggplot2::geom_raster(ggplot2::aes(fill = z)) +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::scale_y_reverse(expand = c(0,0)) +
      ggplot2::scale_fill_distiller(palette = "YlGnBu") +
      ggplot2::labs(x = "Time (years)", y = "Period (years)") +
      ggplot2::guides(fill = "none") +
      ggplot2::geom_line(data = df2, ggplot2::aes(x = x, y = y),
                         linetype = "dashed", color = "red", linewidth = 0.85) +
      ggplot2::stat_contour(ggplot2::aes(z = z), data = df3, breaks = c(-99, 1), color = "black") +
      ggplot2::ggtitle("a)")

    p3 <- ggplot2::ggplot(GWS_gg) +
      ggplot2::theme_light() +
      ggplot2::geom_line(ggplot2::aes(period, GWS)) +
      ggplot2::geom_point(ggplot2::aes(period, GWS), shape = 19, size = 1) +
      ggplot2::geom_line(ggplot2::aes(period, GWS_signif), color = "red", linetype = "dashed", linewidth = 0.85) +
      ggplot2::scale_y_continuous() +
      ggplot2::scale_x_reverse(expand = c(0,0)) +
      ggplot2::coord_flip() +
      ggplot2::labs(y = bquote(Power~(.(variable.unit)^2)), x = "") +
      ggplot2::ggtitle("b)")

    # Save results to file
    p <- p2 + p3 + patchwork::plot_layout(widths = c(2, 1.25))
    ggsave(file.path(output.path, "warm_hist_wavelet_analysis.png"), width=8, height = 6)


  }

  return(list(GWS = GWS,
              GWS_signif = GWS_signif,
              GWS_period = period,
              signif_periods = signif_periods))

}






#
# n1 <- length(variable_org)
# GWS_gg <- tibble::tibble(period = period, GWS = GWS, GWS_signif = GWS_signif)
# var_gg <- tibble::tibble(x = 1:n1, y = variable_org)
#
# df1 <- tibble::as_tibble(t(log2(POWER)), .name_repair = ~ as.character(period)) %>%
#   dplyr::mutate(x = 1:n1) %>%
#   tidyr::gather(key = "y", value = "z", -x) %>%
#   dplyr::mutate_all(as.numeric)
#
# df2 <- tibble::tibble(x = 1:n1, y = coi)
#
# df3 <- tibble::as_tibble(t(sig95), .name_repair = ~ as.character(period)) %>%
#   dplyr::mutate(x = 1:n1) %>%
#   tidyr::gather(key = "y", value = "z", -x) %>%
#   dplyr::mutate_all(as.numeric)
#
# p2 <- ggplot2::ggplot(df1, ggplot2::aes(x = x, y = y)) +
#   ggplot2::theme_light() +
#   ggplot2::geom_raster(ggplot2::aes(fill = z)) +
#   ggplot2::scale_x_continuous(expand = c(0,0)) +
#   ggplot2::scale_y_reverse(expand = c(0,0)) +
#   ggplot2::scale_fill_distiller(palette = "YlGnBu") +
#   ggplot2::labs(x = "Time (years)", y = "Period (years)") +
#   ggplot2::guides(fill = "none") +
#   ggplot2::geom_line(data = df2, ggplot2::aes(x = x, y = y),
#                      linetype = "dashed", color = "red", linewidth = 0.85) +
#   ggplot2::stat_contour(ggplot2::aes(z = z), data = df3, breaks = c(-99, 1), color = "black") +
#   ggplot2::ggtitle("a)")
#
# p3 <- ggplot2::ggplot(GWS_gg) +
#   ggplot2::theme_light() +
#   ggplot2::geom_line(ggplot2::aes(period, GWS)) +
#   ggplot2::geom_point(ggplot2::aes(period, GWS), shape = 19, size = 1) +
#   ggplot2::geom_line(ggplot2::aes(period, GWS_signif), color = "red", linetype = "dashed", linewidth = 0.85) +
#   ggplot2::scale_y_continuous() +
#   ggplot2::scale_x_reverse(expand = c(0,0)) +
#   ggplot2::coord_flip() +
#   ggplot2::labs(y = bquote(Power~(.(variable.unit)^2)), x = "") +
#   ggplot2::ggtitle("b)")
#
# p <- p2 + p3 + patchwork::plot_layout(widths = c(2, 1.25))
# ggplot2::ggsave(file.path(output.path, "wavelet_analysis.png"), plot = p, width = 8, height = 6)
