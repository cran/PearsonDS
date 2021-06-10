\name{pearsonDiagram}
\Rdversion{1.1}
\alias{pearsonDiagram}
\title{
  Regions of Pearson Distribution Types
}
\description{
  2D-Plot of the regions for the different types of Pearson distributions,
  depending on (squared) skewness and kurtosis.
}
\usage{
pearsonDiagram(max.skewness = sqrt(14), max.kurtosis = 24, 
               squared.skewness = TRUE, lwd = 2, legend = TRUE, 
               n = 301)
}
\arguments{
  \item{max.skewness}{
    maximal value for the skewness.
}
  \item{max.kurtosis}{
    maximal value for the kurtosis (not excess kurtosis!).
}
  \item{squared.skewness}{
    plot squared skewness on x-axis (default: \code{TRUE})?
}
  \item{lwd}{
    line width for distributions of type II, III, V, VII.
}
  \item{legend}{
    include legend in the plot (default: \code{TRUE})?
}
  \item{n}{
    number of points for \code{\link{curve}}.
}
}
\value{
  Nothing useful. Function called for its side-effects.
}
\details{
  The label of the x-axis is \eqn{\beta_1}{beta[1]} for squared skewness and
  \eqn{\sqrt{\beta_1}}{sqrt(beta[1])} for regular skewness. The label of the
  y-axis is \eqn{\beta_2}{beta[2]}.
}
\references{
  [1] Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate 
  Distributions, Vol. 1, Wiley Series in Probability and Mathematical Statistics, 
  Wiley
}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de} and
Stefan Klößner \email{S.Kloessner@mx.uni-saarland.de}

}
\seealso{
  \code{\link{PearsonDS-package}}
}
\examples{
## Show me the regions for the different distribution types!
pearsonDiagram()
}
\keyword{ distribution }
