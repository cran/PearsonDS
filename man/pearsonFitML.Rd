\name{pearsonFitML}
\Rdversion{1.1}
\alias{pearsonFitML}
\title{
  Maximum Likelihood Estimation of Pearson Distributions
}
\description{
  This function tries to find the Maximum Likelihood estimator within the 
  Pearson distribution system. ML estimation is done for all sub-classes of
  the distribution system via numerical optimization 
  (with \code{\link{nlminb}}). The sub-class with the
  optimal likelihood function value and the corresponding parameters are
  returned.
}
\usage{
pearsonFitML(x, ...)
}
\arguments{
  \item{x}{
    empirical data (numerical vector) for MLE.
}
  \item{\dots}{
    parameters for \code{\link{nlminb}}.
}
}
\details{
  Starting values for each sub-class are found in a three-step procedure. 
  First, the empirical moments of the input vector are calculated.
  In the second step, the moments are altered, such that the moment restrictions
  for the current sub-class are fulfilled (if necessary), and the method of 
  moments estimator is calculated to obtain starting values for the optimizer.
  In the last step, the starting values are adjusted (if necessary) in order to
  assure that the whole sample lies in the support of the distribution.
}
\note{
  The implementation is VERY preliminary (and slow). No analytical results are 
  used, ie. no analytical solutions for ML estimators and no analytical 
  gradients.
  Most of the distribution types (0, II, III, V, VII) should rather be 
  neglected (for speed reasons), because they will contain the MLE with 
  probability of 0.
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
  \code{\link{pearsonFitM}},
  \code{\link{pearsonMoments}}
}
\examples{
## Generate sample 
DATA <- rpearson(1000,moments=c(mean=1,variance=2,skewness=1,kurtosis=5))
## find Pearson distribution with these parameters
ppar <- pearsonFitML(DATA)
print(unlist(ppar))
## compare with method of moments estimator
print(unlist(pearsonFitM(moments=empMoments(DATA))))
}
\keyword{ distribution }
