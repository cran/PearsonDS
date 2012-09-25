\name{pearsonFitM}
\Rdversion{1.1}
\alias{pearsonFitM}
\title{
  Method of Moments Estimator for Pearson Distributions
}
\description{
  This function calculates the method of moments estimator for Pearson 
  distribution, ie, it generates a Pearson distribution with moments exaktly
  (up to rounding errors) matching the input moments 
  mean, variance, skewness and kurtosis.
}
\usage{
pearsonFitM(mean, variance, skewness, kurtosis, moments)
}
\arguments{
  \item{mean}{
    target mean.
}
  \item{variance}{
    target variance.
}
  \item{skewness}{
    target skewness.
}
  \item{kurtosis}{
    target kurtosis (not excess kurtosis!).
}
  \item{moments}{
    optional vector/list of mean, variance, skewness, kurtosis 
    (not excess kurtosis) in this order. 
    Overrides all other input parameters, if given.
}
}
\value{
    List of parameters for Pearson distribution. First entry gives type 
    of distribution (0 for type 0, 1 for type I, ..., 7 for type VII), remaining
    entries give distribution parameters (depending on distribution type).
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
  \code{\link[PearsonDS]{Pearson}},
  \code{\link{pearsonFitML}},
  \code{\link{pearsonMoments}},
  \code{\link{pearsonMSC}}
}
\examples{
## Define moments of distribution
moments <- c(mean=1,variance=2,skewness=1,kurtosis=5)
## find Pearson distribution with these parameters
ppar <- pearsonFitM(moments=moments)
print(unlist(ppar))
## check moments
pearsonMoments(params=ppar)
}
\keyword{ distribution }
