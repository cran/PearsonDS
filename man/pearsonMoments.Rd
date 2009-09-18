\name{pearsonMoments}
\Rdversion{1.1}
\alias{pearsonMoments}
\title{
  Moments of Pearson Distribution
}
\description{
  Calculates the first four moments (mean, variance, skewness, kurtosis)
  of a Pearson distribution.
}
\usage{
pearsonMoments(params, moments)
}
\arguments{
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
}
\value{
  Named vector of length 4 containing \code{mean}, \code{variance}, 
  \code{skewness} and \code{kurtosis} (in this order).
}
\note{
  Optional parameter \code{moments} is merely for testing purposes. Of course,
  \code{pearsonMoments} should reproduce its input (when neglecting rounding
  errors) if \code{moments} is given.
  
  \code{kurtosis} is the kurtosis of the distribution, not the excess kurtosis!
}
\references{
  [1] Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate 
  Distributions, Vol. 1, Wiley Series in Probability and Mathematical Statistics, 
  Wiley

  [2] Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate 
  Distributions, Vol. 2, Wiley Series in Probability and Mathematical Statistics, 
  Wiley
}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\seealso{
  \code{\link{PearsonDS-package}},
  \code{\link{Pearson0}}, \code{\link{PearsonI}}, \code{\link{PearsonII}}, 
  \code{\link{PearsonIII}}, \code{\link{PearsonIV}}, \code{\link{PearsonV}}, 
  \code{\link{PearsonVI}}, \code{\link{PearsonVII}}, 
}
\examples{
## Define moments of distribution
moments <- c(mean=1,variance=2,skewness=1,kurtosis=5)
## Are the moments reproduced?
pearsonMoments(moments=moments)
}
\keyword{ distribution }
