% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateData.R
\name{initializeGphmm}
\alias{initializeGphmm}
\title{Initial parameters for gphmm.}
\usage{
initializeGphmm()
}
\description{
\code{initializeGphmm} returns initial set of parameters for gphmm.
}
\examples{
param <- initializeGphmm ()
param[['qR']] <- c(0.1, 0.3, 0.3, 0.1)
}
