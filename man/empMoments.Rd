\name{empMoments}
\Rdversion{1.1}
\alias{empMoments}
\title{
  Empirical Moments
}
\description{
  Calculates the first four empirical moments (mean, variance, skewness, 
  kurtosis) of a numeric vector.
}
\usage{
empMoments(x)
}
\arguments{
  \item{x}{
    (numeric) vector containing the data set.
}
}
\value{
  Named vector of length 4 containing \code{mean}, \code{variance}, 
  \code{skewness} and \code{kurtosis} (in this order).
}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\seealso{
  \code{\link{PearsonDS-package}}
}
\examples{
## Generate sample with given (theoretical) moments
DATA <- rpearson(25000,moments=c(mean=1,variance=2,skewness=1,kurtosis=5))
## Calculate corresponding empirical moments
empMoments(DATA)
}
\keyword{ distribution }
