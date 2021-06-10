\name{PearsonI}
\Rdversion{1.1}
\alias{PearsonI}
\alias{dpearsonI}
\alias{ppearsonI}
\alias{qpearsonI}
\alias{rpearsonI}
\title{
  The Pearson Type I (aka Beta) Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type I (aka Beta) distribution.
}
\usage{

dpearsonI(x, a, b, location, scale, params, log = FALSE)

ppearsonI(q, a, b, location, scale, params, lower.tail = TRUE, 
          log.p = FALSE)

qpearsonI(p, a, b, location, scale, params, lower.tail = TRUE, 
          log.p = FALSE)

rpearsonI(n, a, b, location, scale, params)
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
    first shape parameter of Pearson type I distribution.
}
  \item{b}{
    second shape parameter of Pearson type I distribution.
}
  \item{location}{
    location parameter of Pearson type I distribution.
}
  \item{scale}{
    scale parameter of Pearson type I distribution.
}
  \item{params}{
    vector/list of length 4 containing parameters \code{a}, \code{b},
    \code{location}, \code{scale} for Pearson type I distribution 
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
  Essentially, Pearson type I distributions are (location-scale transformations 
  of) Beta distributions, the above 
  functions are thus simple wrappers for \code{dbeta}, \code{pbeta}, 
  \code{qbeta} and \code{rbeta} contained in package \code{stats}.
  The probability density function with parameters \code{a}, \code{b},  
  \code{scale}\eqn{=s} and \code{location}\eqn{=\lambda}{=lambda} 
  is given by
  \deqn{f(x)=\frac{1}{|s|}\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}\left(\frac{x-\lambda}{s}
  \right)^{a-1}\left(1-\frac{x-\lambda}{s}\right)^{b-1}}{f(x)=1/|s|*Gamma(a+b)/
  (Gamma(a)Gamma(b))
  ((x-lambda)/s)^(a-1)(1-((x-lambda)/s))^(b-1)}
  for \eqn{a>0}, \eqn{b>0}, \eqn{s\ne 0}{s<>0}, 
  \eqn{0<\frac{x-\lambda}{s}<1}{0<(x-lambda)/s<1}.
}
\value{
  \code{dpearsonI} gives the density, \code{ppearsonI} gives the 
  distribution function, \code{qpearsonI} gives the quantile function, 
  and \code{rpearsonI} generates random deviates.
}
\references{
  See the references in \code{\link{Beta}}.
}
\author{
  Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\note{
  Negative values for \code{scale} are not excluded, albeit negative skewness 
  is usually obtained by switching \code{a} and \code{b} 
  (such that \code{a}>\code{b}) and not by using negative values for 
  \code{scale} (and \code{a}<\code{b}).
}
\seealso{
  \code{\link{Beta}},
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type I parameter set with a=2, b=3, location=1, scale=2
pIpars <- list(a=2, b=3, location=1, scale=2)
## calculate probability density function
dpearsonI(seq(1,3,by=0.5),params=pIpars)
## calculate cumulative distribution function
ppearsonI(seq(1,3,by=0.5),params=pIpars)
## calculate quantile function
qpearsonI(seq(0.1,0.9,by=0.2),params=pIpars)
## generate random numbers
rpearsonI(5,params=pIpars)
}
\keyword{ distribution }
