\name{matchMoments}
\alias{matchMoments}
\title{
  Fitting (Incomplete) Set of Moments to Given Distribution Type
}
\description{
  For a given incomplete (skewness or kurtosis are missing) set of moments,
  the complete set of moments (mean, variance, skewness, and kurtosis) is
  calculated, using a given distribution type (if possible). Either the
  complete set of moments or the distribution parameters are returned.
}
\usage{
matchMoments(mean, variance, skewness = NA, kurtosis = NA, type, moments,
             skewness.sign = c("+", "-"), return.distribution = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mean}{
    target mean.
}
  \item{variance}{
    target variance.
}
  \item{skewness}{
    target skewness (maybe \code{NA}, which is the default).
}
  \item{kurtosis}{
    target kurtosis (not excess kurtosis, maybe \code{NA}, which is the default).
}
  \item{type}{
    required distribution type (either one of the numbers 0, 1, ..., 7 or one of
    \code{"0"}, \code{"I"}, \code{"II"}, \code{"III"}, \code{"IV"}, \code{"V"},
    \code{"VI"}, or  \code{"VII"}).
}
  \item{moments}{
    optional vector/list of mean, variance, skewness, kurtosis (not excess
    kurtosis) in this order. Overrides the input parameters \code{mean},
    \code{variance}, \code{skewness}, \code{kurtosis} (if given).
}
  \item{skewness.sign}{
    \code{"+""} or \code{"-""}. Gives the sign of the resulting skewness if
    the input parameter \code{skewness} is \code{NA} and the resulting skewness
    is not zero (default: \code{"+"}).
}
  \item{return.distribution}{
    \code{TRUE} or \code{FALSE}. If set to \code{FALSE} (default), the
    complete set of moments is returned, if set to \code{TRUE}, the list of
    parameters for the corresponging Pearson distribution is returned.
}
}
\value{
  If \code{return.distribution==TRUE}: list of parameters for Pearson
  distribution. First entry gives type of distribution (0 for type 0,
  1 for type I, ..., 7 for type VII), remaining entries give distribution
  parameters (depending on distribution type).

  If \code{return.distribution==FALSE}: numeric vector with named elements
  \code{mean}, \code{variance}, \code{skewness}, \code{kurtosis} corresponding
  to a Pearson distribution of type \code{type}.
}

\note{
  If the specified subset of moments does not match the distribution type or
  if too many moments are missing, an error is issued.
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
matchMoments(mean=0,variance=1,kurtosis=4.5,type=3,return.distribution = TRUE)
matchMoments(mean=0,variance=1,kurtosis=4.5,type="III")
matchMoments(mean=0,variance=1,kurtosis=4.5,type="III",skewness.sign="-")

matchMoments(mean=0,variance=1,skewness=-2,type="III",return.distribution = TRUE)
pearsonFitM(moments=matchMoments(mean=0,variance=1,skewness=-2,type="III"))
matchMoments(mean=0,variance=1,skewness=-2,type="III")
}

\keyword{ distribution }
