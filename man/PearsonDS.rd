\name{PearsonDS-package}
\Rdversion{1.1}
\alias{PearsonDS-package}
\alias{PearsonDS}
\docType{package}
\title{
Pearson Distribution System
}
\description{
  Implementation of the \code{d,p,q,r} function family, calculation of moments, 
  and fitting via (empirical) moment matching as well as maximum likelihood 
  method for the Pearson distribution system.
}
\section{Warning}{  
  If at all possible, package \code{gsl} should be installed. 
  In this case, the functions for Pearson type IV 
  distributions make use of \code{lngamma_complex} 
  (see \code{\link[gsl]{Gamma}}).
  If package \code{gsl} is not installed, some
  calculations for Pearson type IV distributions with (more or less) extreme 
  parameters (ie, big \code{nu} and/or \code{m}) may slow down by factors of
  more than 1000.
}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de} 
and Stefan Klößner \email{S.Kloessner@mx.uni-saarland.de}

Maintainer: Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\references{
  [1] Abramowitz, M. and I. A. Stegun (1972) 
  \emph{Handbook of mathematical functions}, 
  National Bureau of Standards, Applied Mathematics Series - 55, Tenth Printing,
  Washington D.C.

  [2] Heinrich, J. (2004) \emph{A Guide to the Pearson Type IV Distribution}, 
  Univ. Pennsylvania, Philadelphia, Tech. Rep. CDF/Memo/Statistics/Public/6820
  \url{http://www-cdf.fnal.gov/physics/statistics/notes/cdf6820_pearson4.pdf}

  [3] Hida, Y., X. S. Li and D. H. Bailey (2000) 
  \emph{Algorithms for quad-double precision floating point arithmetic}, 
  Lawrence Berkeley National Laboratory. Paper LBNL-48597.

  [4] Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate 
  Distributions, Vol. 1, Wiley Series in Probability and Mathematical Statistics, 
  Wiley

  [5] Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate 
  Distributions, Vol. 2, Wiley Series in Probability and Mathematical Statistics, 
  Wiley

  [6] Willink, R. (2008) \emph{A Closed-form Expression for the Pearson Type IV 
  Distribution Function}, Aust. N. Z. J. Stat. 50 (2), pp. 199-205
}
\keyword{ package }
\keyword{ distribution } 
\seealso{
  \code{\link[PearsonDS]{Pearson}} for \code{d,p,q,r} function family for Pearson 
  distributions, 
  \code{\link{pearsonFitM}} and \code{\link{pearsonFitML}} for fitting 
  Pearson distributions, 
  \code{\link{pearsonMSC}} for model selection,
  \code{\link{pearsonMoments}} for calculation of (first four) moments.
}
\examples{
## see documentation of individual functions
}
