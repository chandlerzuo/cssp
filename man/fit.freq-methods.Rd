% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/post_model.R
\docType{methods}
\name{fit.freq}
\alias{fit.freq}
\alias{fit.freq,CSSPFit-method}
\title{Compute the estimated frequency for ChIP counts based on the CSSP model.}
\usage{
fit.freq(fit, chip)

\S4method{fit.freq}{CSSPFit}(fit, chip)
}
\arguments{
\item{fit}{A \link{CSSPFit-class} object for the fitted CSSP model.}

\item{chip}{A \link{numeric} vector of ChIP sample bin counts.}
}
\value{
A \link{data.frame} object containing
\tabular{l}{
count The counts of each bin.\\
freq The ChIP data frequency at this count value.\\
freq.est  The estimated frequency using the posterior distributions of the bin-level poisson intensities.}
}
\description{
Compute the estimated frequency for ChIP counts based on the CSSP model.
}
\examples{
data( sampleFit )
data( bin.data )
fit.freq( sampleFit, chip = bin.data@tagCount )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

