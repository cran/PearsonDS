\name{PearsonVI}
\Rdversion{1.1}
\alias{PearsonVI}
\alias{dpearsonVI}
\alias{ppearsonVI}
\alias{qpearsonVI}
\alias{rpearsonVI}
\title{
  The Pearson Type VI (aka Beta Prime) Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type VI (aka Beta prime) distribution.
}
\usage{

dpearsonVI(x, a, b, location, scale, params, log = FALSE)

ppearsonVI(q, a, b, location, scale, params, lower.tail = TRUE, 
           log.p = FALSE)

qpearsonVI(p, a, b, location, scale, params, lower.tail = TRUE, 
           log.p = FALSE)

rpearsonVI(n, a, b, location, scale, params)
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
    first shape parameter of Pearson type VI distribution.
}
  \item{b}{
    second shape parameter of Pearson type VI distribution.
}
  \item{location}{
    location parameter of Pearson type VI distribution.
}
  \item{scale}{
    scale parameter of Pearson type VI distribution.
}
  \item{params}{
    vector/list of length 4 containing parameters \code{a}, \code{b},
    \code{location}, \code{scale} for Pearson type VI distribution 
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
  Pearson type VI distributions are (location-scale transformations 
  of) Beta prime distributions, and Beta prime distributions are 
  scaled F-distributions. The above functions are 
  thus implemented via calls to \code{df}, \code{pf}, 
  \code{qf} and \code{rf} (contained in package \code{stats}).
  The probability density function with parameters \code{a}, \code{b},  
  \code{scale}\eqn{=s} and \code{location}\eqn{=\lambda}{=lambda} 
  is given by
  \deqn{f(x)=\frac{\Gamma(a+b)}{|s|\Gamma(a)\Gamma(b)}\left(\frac{x-\lambda}{s}
  \right)^{a-1}\left(1+\frac{x-\lambda}{s}\right)^{-a-b}}{f(x)=Gamma(a+b)/
  (|s|Gamma(a)Gamma(b))((x-lambda)/s)^(a-1)(1+(x-lambda)/s)^(-a-b)}
  for \eqn{a>0}, \eqn{b>0}, \eqn{s\ne 0}{s<>0}, 
  \eqn{\frac{x-\lambda}{s}>0}{(x-lambda)/s>0}.
}
\value{
  \code{dpearsonVI} gives the density, \code{ppearsonVI} gives the 
  distribution function, \code{qpearsonVI} gives the quantile function, 
  and \code{rpearsonVI} generates random deviates.
}
\references{
  See the references in \code{\link{FDist}}.
}
\author{
  Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\note{
  Negative values for \code{scale} are permitted to allow for negative skewness.
}
\seealso{
  \code{\link{FDist}},
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type VI parameter set with a=2, b=3, location=1, scale=2
pVIpars <- list(a=2, b=3, location=1, scale=2)
## calculate probability density function
dpearsonVI(seq(1,6,by=1),params=pVIpars)
## calculate cumulative distribution function
ppearsonVI(seq(1,6,by=1),params=pVIpars)
## calculate quantile function
qpearsonVI(seq(0.1,0.9,by=0.2),params=pVIpars)
## generate random numbers
rpearsonVI(5,params=pVIpars)
}
\keyword{ distribution }
