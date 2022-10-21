\name{Pearson0}
\Rdversion{1.1}
\alias{Pearson0}
\alias{dpearson0}
\alias{ppearson0}
\alias{qpearson0}
\alias{rpearson0}
\title{
  The Pearson Type 0 (aka Normal) Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type 0 (aka Normal) distribution.
}
\usage{

dpearson0(x, mean, sd, params, log = FALSE)

ppearson0(q, mean, sd, params, lower.tail = TRUE, log.p = FALSE)

qpearson0(p, mean, sd, params, lower.tail = TRUE, log.p = FALSE)

rpearson0(n, mean, sd, params)
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
  \item{mean}{
    location parameter (and expectation)
}
  \item{sd}{
    scale parameter (and standard deviation)
}
  \item{params}{
    optional vector/list containing distribution parameters \code{mean} and 
    \code{sd} (in this order!). Overrides parameters \code{mean} and \code{sd},
    if given.
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
  Distributions of type 0 have been added to the Pearson Distribution system 
  in order to have the normal distributions not only nested as limits of other
  distribution types. 
  The functions are only wrappers for \code{dnorm}, \code{pnorm}, \code{qnorm} 
  and \code{rnorm} contained in package \code{stats}.
}
\value{
  \code{dpearson0} gives the density, \code{ppearson0} gives the distribution 
  function, \code{qpearson0} gives the quantile function, and \code{rpearson0} 
  generates random deviates.
}
\references{
  See the references in \code{\link{Normal}}.
}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\seealso{
  \code{\link{Normal}},
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type 0 parameter set with mean=-1, sd=2
p0pars <- list(mean=-1, sd=2)
## calculate probability density function
dpearson0(-4:1,params=p0pars)
## calculate cumulative distribution function
ppearson0(-4:1,params=p0pars)
## calculate quantile function
qpearson0(seq(0.1,0.9,by=0.2),params=p0pars)
## generate random numbers
rpearson0(5,params=p0pars)
}

\keyword{ distribution }
