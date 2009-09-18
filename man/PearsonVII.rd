\name{PearsonVII}
\Rdversion{1.1}
\alias{PearsonVII}
\alias{dpearsonVII}
\alias{ppearsonVII}
\alias{qpearsonVII}
\alias{rpearsonVII}
\title{
  The Pearson Type VII (aka Student's t) Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type VII (aka Student's t) distribution.
}
\usage{

dpearsonVII(x, df, location, scale, params, log = FALSE)

ppearsonVII(q, df, location, scale, params, lower.tail = TRUE, 
            log.p = FALSE)

qpearsonVII(p, df, location, scale, params, lower.tail = TRUE, 
            log.p = FALSE)

rpearsonVII(n, df, location, scale, params)
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
  \item{df}{
    degrees of freedom of Pearson type VII distribution
}
  \item{location}{
    location parameter of Pearson type VII distribution.
}
  \item{scale}{
    scale parameter of Pearson type VII distribution.
}
  \item{params}{
    vector/list of length 3 containing parameters \code{df}, 
    \code{location}, \code{scale} for Pearson type VII distribution 
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
  The Pearson type VII distribution is a simple (location-scale) transformation
  of the well-known Student's t distribution; the probability density function
  with parameters \code{df}\eqn{=n}, \code{location}\eqn{=\lambda}{=lambda} and
  \code{scale}\eqn{=s} is given by
  \deqn{f(x) = \frac{\Gamma(\frac{n+1}{2})}{\sqrt{n \pi} \Gamma(\frac{n}{2})} 
               \left(1 + \frac{(\frac{x-\lambda}{s})^2}{n}\right)^{-\frac{n+1}{2}}}{f(x) 
               = Gamma((n+1)/2) / (sqrt(n pi) Gamma(n/2)) 
               (1 + ((x-lambda)/s)^2/n)^-((n+1)/2)}
  for \eqn{s\ne 0}{s<>0}.
  The above functions are thus only wrappers for \code{dt}, 
  \code{pt}, \code{qt} and \code{rt} contained in package 
  \code{stats}.
}
\value{
  \code{dpearsonVII} gives the density, \code{ppearsonVII} gives the 
  distribution function, \code{qpearsonVII} gives the quantile function, 
  and \code{rpearsonVII} generates random deviates.
}
\references{
  See the references in \code{\link{TDist}}.
}
\author{
  Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\seealso{
  \code{\link{TDist}},
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type VII parameter set with df=7, location=1, scale=1
pVIIpars <- list(df=7, location=1, scale=1)
## calculate probability density function
dpearsonVII(-2:4,params=pVIIpars)
## calculate cumulative distribution function
ppearsonVII(-2:4,params=pVIIpars)
## calculate quantile function
qpearsonVII(seq(0.1,0.9,by=0.2),params=pVIIpars)
## generate random numbers
rpearsonVII(5,params=pVIIpars)
}
\keyword{ distribution }
