% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waveletARIMA.R
\name{waveletARIMA}
\alias{waveletARIMA}
\title{Wavelet Autoregressive Modeling (WARM)}
\usage{
waveletARIMA(
  wavelet.components = NULL,
  sim.year.num = NULL,
  sim.num = 1000,
  seed = NULL
)
}
\arguments{
\item{wavelet.components}{A list object, with different components corresponding to low-frequency signals and the noise}

\item{sim.year.num}{A numeric value defining the desired length of simulated annual time-series}

\item{sim.num}{A numeric value defining the number of synthetic series to be produced.}

\item{seed}{random seed value}
}
\description{
Wavelet Autoregressive Modeling (WARM)
}
