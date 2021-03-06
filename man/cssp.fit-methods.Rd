% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/fit_model.R
\docType{methods}
\name{cssp.fit}
\alias{cssp.fit}
\alias{cssp.fit,BinData-method}
\alias{cssp.fit,data.frame-method}
\alias{cssp.fit,data.frame-method,BinData-method}
\title{Fit the CSSP Model.}
\usage{
cssp.fit(dat, method = "mde", p1 = 0.5, p2 = 0.99, beta.init = NULL,
  e0.init = 0.9, e0.lb = 0.5, ngc = 9, nite = 50, tol = 0.01,
  useGrid = FALSE, nsize = NULL, ncomp = 2, nonpa = FALSE,
  zeroinfl = FALSE, seed = NULL)

\S4method{cssp.fit}{data.frame}(dat, method = "mde", p1 = 0.5, p2 = 0.99,
  beta.init = NULL, e0.init = 0.9, e0.lb = 0.5, ngc = 9, nite = 50,
  tol = 0.01, useGrid = FALSE, nsize = NULL, ncomp = 2, nonpa = FALSE,
  zeroinfl = FALSE, seed = NULL)

\S4method{cssp.fit}{BinData}(dat, method = "mde", p1 = 0.5, p2 = 0.99,
  beta.init = NULL, e0.init = 0.9, e0.lb = 0.5, ngc = 9, nite = 50,
  tol = 0.01, useGrid = FALSE, nsize = NULL, ncomp = 2, nonpa = FALSE,
  zeroinfl = FALSE, seed = NULL)
}
\arguments{
\item{dat}{A \link{data.frame} or \link{BinData-class} object containing bin-level chip, input, M and GC information. For the data.frame object, the columns must contain "chip", "input", "M". For BinData object, the slots must contain "tagCount", "input", "M". If "GC" is not provided, model will be fitted without using gc-Content scores.}

\item{method}{A \link{character} indicating the method of fitting algorithm to be used. "mde" (Default) - minimum distance estimation; "gem" - the generalized EM method.}

\item{p1}{The \link{numeric} value for the lower bound for the p-value region where the p-values are assumed to be uniformly distributed. Default: 0.5.}

\item{p2}{The \link{numeric} value for the upper bound for the p-value region where the p-values are assumed to be uniformly distributed. Default: 0.99.}

\item{beta.init}{The \link{numeric} value for the initializing the size parameter for the background model of the ChIP sample. If "NULL", the size parameter of the fitted input sample model is used.}

\item{e0.init}{The \link{numeric} value for initializing parameter e0. Default: 0.9.}

\item{e0.lb}{The \link{numeric} value for the lower bound of parameter e0. Default is 0.5. This parameter is recommended to be set according to the p-value plot.}

\item{ngc}{An \link{integer} value for the number of knots used in the spline model for the gc covariate. Default: 9.}

\item{nite}{An \link{integer} value for the maximum number of iterations taken. Default: 50.}

\item{tol}{A \link{numeric} value for the tolerance for convergence. Default: 1e-3.}

\item{useGrid}{A \link{logical} value indicating whether the gridding method is used. If TRUE, the covariate space is grided adaptively. This trims down the sample size for fitting the regression model when the data contains too many observations, and is suggested for genome-wide analysis. Default: FALSE.}

\item{nsize}{A \link{numeric} value for the number of bins to be randomly chosen in estimating the normalizatiing parameters. If Null (default), all bins are used in normalization. For genome wide analysis, nsize=5000 is suggested.}

\item{ncomp}{A \link{numeric} value for the number of signal components.}

\item{nonpa}{A \link{logical} value indicating whether a nonparametric model for the background ChIP sample and the input sample is fitted.}

\item{zeroinfl}{A \link{logical} value indicating whether a zero-inflated negative binomial model is fitted for the ChIP background.}

\item{seed}{A \link{numeric} value for the seed of generating random variables. Default: NULL. Users should specify this value for generating exactly reproducible results.}
}
\value{
\link{CSSPFit-class} A CSSPFit object.
}
\description{
Fit the CSSP Model.
}
\details{
The current version of cssp.fit has implemented the following method.\cr
The "method" argument specifies the method to estimate the normalization models for the ChIP background from the input data. "mde" uses minimum distance estimation, "gem" uses generalized E-M estimation.\cr
The 'nonpa' argument specifies whether a glm model is used. If "nonpa" is FALSE, a GLM is used to fit the input data. If "nonpa" is TRUE, the mean response within each grid is taken as the predict. These two arguments enables the analysis for genome-wide data. In this case, "nsize" grids are used.\cr
If "nonpa" is FALSE, then "useGrid" specifies whether the covariate space is grided adaptively, and the mean values within each grid is used for regression.\cr
If "nonpa" is TRUE, "zeroinfl" specifies whether a zero-inflation model for the background is used. This is useful for low-depth ChIP data, where too many bins have zero count.
}
\examples{
data( bin.data )
cssp.fit( bin.data )
cssp.fit( bin.data, method = "gem" )
data( bindata.chr1 )
cssp.fit( bindata.chr1 )
cssp.fit( bindata.chr1, method = "gem", ngc = 1 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

