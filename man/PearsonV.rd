\name{PearsonV}
\Rdversion{1.1}
\alias{PearsonV}
\alias{dpearsonV}
\alias{ppearsonV}
\alias{qpearsonV}
\alias{rpearsonV}
\title{
  The Pearson Type V (aka Inverse Gamma) Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type V (aka Inverse Gamma) distribution.
}
\usage{

dpearsonV(x, shape, location, scale, params, log = FALSE)

ppearsonV(q, shape, location, scale, params, lower.tail = TRUE, 
          log.p = FALSE)

qpearsonV(p, shape, location, scale, params, lower.tail = TRUE, 
          log.p = FALSE)

rpearsonV(n, shape, location, scale, params)
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
  \item{shape}{
    shape parameter of Pearson type V distribution.
}
  \item{location}{
    location parameter of Pearson type V distribution.
}
  \item{scale}{
    scale parameter of Pearson type V distribution.
}
  \item{params}{
    vector/list of length 3 containing parameters \code{shape}, 
    \code{location}, \code{scale} for Pearson type V distribution 
    (in this order!). 
}
  \item{log, log.p}{
    logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
}
  \item{lower.tail}{
    logical; if \code{TRUE}, probabilities are \eqn{P[X\le x]}{P[X<=x]}, 
    otherwise, \eqn{P[X>x]}{P[X>x]}.
}
}
\details{
  The Pearson type V distributions are essentially Inverse Gamma distributions.
  Thus, all functions are implemented via calls to the corresponding functions 
  for Gamma distributions, ie. \code{dgamma}, \code{pgamma}, \code{qgamma} 
  and \code{rgamma} in package \code{stats}.
  Negative \code{scale} parameters 
  (which reflect the distribution at \code{location}) are 
  permitted to allow for negative skewness.
  The probability density function with parameters \code{shape}\eqn{=a},  
  \code{scale}\eqn{=s} and \code{location}\eqn{=\lambda}{=lambda} 
  is given by
  \deqn{f(x)= \frac{|s|^a }{\Gamma(a)} |x-\lambda|^{-a-1} 
              e^{-\frac{s}{x-\lambda}}}{f(x)= |s|^a/Gamma(a)
              |x-lambda|^(-a-1) e^-(s/(x-lambda))}
  for \eqn{s\ne 0}{s<>0}, \eqn{a>0} and 
  \eqn{\frac{s}{x-\lambda}> 0}{s/(x-lambda)>0}.
}
\value{
  \code{dpearsonV} gives the density, \code{ppearsonV} gives the 
  distribution function, \code{qpearsonV} gives the quantile function, 
  and \code{rpearsonV} generates random deviates.
}
\references{
  See the references in \code{\link{GammaDist}}.
}
\author{
  Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\note{
  Since package version 0.98, the parameter \code{scale} corresponds to the 
  usual scale parameter of the Inverse Gamma distribution (not the reciprocal 
  value, which was implemented [incorrectly!] until package version 0.97). 
}
\seealso{
  \code{\link{GammaDist}},
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type V parameter set with shape=3, location=1, scale=-2
pVpars <- list(shape=3, location=1, scale=-2)
## calculate probability density function
dpearsonV(-4:1,params=pVpars)
## calculate cumulative distribution function
ppearsonV(-4:1,params=pVpars)
## calculate quantile function
qpearsonV(seq(0.1,0.9,by=0.2),params=pVpars)
## generate random numbers
rpearsonV(5,params=pVpars)
}
\keyword{ distribution }
