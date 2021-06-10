\name{PearsonIII}
\Rdversion{1.1}
\alias{PearsonIII}
\alias{dpearsonIII}
\alias{ppearsonIII}
\alias{qpearsonIII}
\alias{rpearsonIII}
\title{
  The Pearson Type III (aka Gamma) Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type III (aka Gamma) distribution.
}
\usage{

dpearsonIII(x, shape, location, scale, params, log = FALSE)

ppearsonIII(q, shape, location, scale, params, lower.tail = TRUE, 
            log.p = FALSE)

qpearsonIII(p, shape, location, scale, params, lower.tail = TRUE, 
            log.p = FALSE)

rpearsonIII(n, shape, location, scale, params)
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
    shape parameter of Pearson type III distribution.
}
  \item{location}{
    location parameter of Pearson type III distribution.
}
  \item{scale}{
    scale parameter of Pearson type III distribution.
}
  \item{params}{
    vector/list of length 3 containing parameters \code{shape}, 
    \code{location}, \code{scale} for Pearson type III distribution 
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
  Essentially, the above functions are wrappers for \code{dgamma}, 
  \code{pgamma}, \code{qgamma} and \code{rgamma} contained in package 
  \code{stats}.
  As a minor (but important) extension, negative \code{scale} parameters 
  (which reflect the distribution at \code{location}) are 
  permitted to allow for negative skewness.
  The probability density function with parameters \code{shape}\eqn{=a},  
  \code{scale}\eqn{=s} and \code{location}\eqn{=\lambda}{=lambda} 
  is thus given by
  \deqn{f(x)= \frac{1}{|s|^a \Gamma(a)} |x-\lambda|^{a-1} 
              e^{-\frac{x-\lambda}{s}}}{f(x)= 1/(|s|^a Gamma(a)) 
              |x-lambda|^(a-1) e^-((x-lambda)/s)}
  for \eqn{s\ne 0}{s<>0}, \eqn{a>0} and 
  \eqn{\frac{x-\lambda}{s}\ge 0}{(x-lambda)/s>=0}.
}
\value{
  \code{dpearsonIII} gives the density, \code{ppearsonIII} gives the 
  distribution function, \code{qpearsonIII} gives the quantile function, 
  and \code{rpearsonIII} generates random deviates.
}
\references{
  See the references in \code{\link{GammaDist}}.
}
\author{
  Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\seealso{
  \code{\link{GammaDist}},
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type III parameter set with shape=3, location=1, scale=-2
pIIIpars <- list(shape=3, location=1, scale=-0.5)
## calculate probability density function
dpearsonIII(-4:1,params=pIIIpars)
## calculate cumulative distribution function
ppearsonIII(-4:1,params=pIIIpars)
## calculate quantile function
qpearsonIII(seq(0.1,0.9,by=0.2),params=pIIIpars)
## generate random numbers
rpearsonIII(5,params=pIIIpars)
}
\keyword{ distribution }
