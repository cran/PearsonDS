\name{PearsonII}
\Rdversion{1.1}
\alias{PearsonII}
\alias{dpearsonII}
\alias{ppearsonII}
\alias{qpearsonII}
\alias{rpearsonII}
\title{
  The Pearson Type II (aka Symmetric Beta) Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type II (aka symmetric Beta) distribution.
}
\usage{

dpearsonII(x, a, location, scale, params, log = FALSE)

ppearsonII(q, a, location, scale, params, lower.tail = TRUE, 
           log.p = FALSE)

qpearsonII(p, a, location, scale, params, lower.tail = TRUE, 
           log.p = FALSE)

rpearsonII(n, a, location, scale, params)
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
  \item{a}{
    shape parameter of Pearson type II distribution.
}
  \item{location}{
    location parameter of Pearson type II distribution.
}
  \item{scale}{
    scale parameter of Pearson type II distribution.
}
  \item{params}{
    vector/list of length 3 containing parameters \code{a},
    \code{location}, \code{scale} for Pearson type II distribution 
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
  Essentially, Pearson type II distributions are (location-scale transformations 
  of) symmetric Beta distributions, the above 
  functions are thus simple wrappers for \code{dbeta}, \code{pbeta}, 
  \code{qbeta} and \code{rbeta} contained in package \code{stats}.
  The probability density function with parameters \code{a},  
  \code{scale}\eqn{=s} and \code{location}\eqn{=\lambda}{=lambda} 
  is given by
  \deqn{f(x)=\frac{\Gamma(2a)}{\Gamma(a)^2}\left(\frac{x-\lambda}{s}\cdot
        \left(1-\frac{x-\lambda}{s}\right)\right)^{a-1}}{f(x)=Gamma(2a)/
        (Gamma(a)^2)
        (((x-lambda)/s)*(1-((x-lambda)/s)))^(a-1)}
  for \eqn{a>0}, \eqn{s\ne 0}{s<>0}, 
  \eqn{0<\frac{x-\lambda}{s}<1}{0<(x-lambda)/s<1}.
}
\value{
  \code{dpearsonII} gives the density, \code{ppearsonII} gives the 
  distribution function, \code{qpearsonII} gives the quantile function, 
  and \code{rpearsonII} generates random deviates.
}
\references{
  See the references in \code{\link{Beta}}.
}
\author{
  Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\seealso{
  \code{\link{Beta}},
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type II parameter set with a=2, location=1, scale=2
pIIpars <- list(a=2, location=1, scale=2)
## calculate probability density function
dpearsonII(seq(1,3,by=0.5),params=pIIpars)
## calculate cumulative distribution function
ppearsonII(seq(1,3,by=0.5),params=pIIpars)
## calculate quantile function
qpearsonII(seq(0.1,0.9,by=0.2),params=pIIpars)
## generate random numbers
rpearsonII(5,params=pIIpars)
}
\keyword{ distribution }
