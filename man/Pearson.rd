\name{Pearson}
\Rdversion{1.1}
\alias{dpearson}
\alias{ppearson}
\alias{qpearson}
\alias{rpearson}
\title{
  The Pearson Distribution System
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson distribution system.
}
\usage{

dpearson(x, params, moments, log = FALSE, ...)

ppearson(q, params, moments, lower.tail = TRUE, log.p = FALSE, ...)

qpearson(p, params, moments, lower.tail = TRUE, log.p = FALSE, ...)

rpearson(n, params, moments, ...)
}
\arguments{

  \item{x, q}{
    vector of quantiles.
}
  \item{p}{
    vector of probabilities.
}
  \item{n}{
    number of observations.
}
  \item{params}{
    vector/list of parameters for Pearson distribution. First entry gives type 
    of distribution (0 for type 0, 1 for type I, ..., 7 for type VII), remaining
    entries give distribution parameters (depending on distribution type).
}
  \item{moments}{
    optional vector/list of mean, variance, skewness, kurtosis 
    (not excess kurtosis). 
    Overrides \code{params} with corresponding pearson distribution, if given.
}
  \item{log, log.p}{
    logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
}
  \item{lower.tail}{
    logical; if \code{TRUE}, probabilities are \eqn{P[X\le x]}{P[X<=x]}, 
    otherwise, \eqn{P[X>x]}{P[X>x]}.
}
  \item{\dots}{
    further parameters for underlying functions (currently only used for
    distributions of type IV).
}
}
\details{
  These are the wrapper functions for the (d,p,q,r)-functions of the Pearson
  distribution system sub-classes.
}
\value{
  \code{dpearson} gives the density, \code{ppearson} gives the distribution 
  function, \code{qpearson} gives the quantile function, and \code{rpearson} 
  generates random deviates.
}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\seealso{
  \code{\link{PearsonDS-package}},
  \code{\link{Pearson0}}, \code{\link{PearsonI}}, \code{\link{PearsonII}}, 
  \code{\link{PearsonIII}}, \code{\link{PearsonIV}}, \code{\link{PearsonV}}, 
  \code{\link{PearsonVI}}, \code{\link{PearsonVII}}, \code{\link{pearsonFitM}},
  \code{\link{pearsonFitML}} 
}
\examples{
## Define moments of distribution
moments <- c(mean=1,variance=2,skewness=1,kurtosis=5)
## Generate some random variates
rpearson(5,moments=moments)
## evaluate distribution function
ppearson(seq(-2,3,by=1),moments=moments)
## evaluate density function
dpearson(seq(-2,3,by=1),moments=moments)
## evaluate quantile function
qpearson(seq(0.1,0.9,by=0.2),moments=moments)
}
\keyword{ distribution }
