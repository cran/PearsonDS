\name{pearsonMSC}
\alias{pearsonMSC}
\title{
  Log-Likelihoods and Model Selection Criteria for different Pearson 
  distribution types
}
\description{
  This function performs (as \code{\link{pearsonFitML}}) an ML estimation 
  for all sub-classes of the Pearson distribution system via numerical 
  optimization (with \code{\link{nlminb}}) for model selection purposes.
  Apart from calculating the log-likelihood values as well as the values of 
  some common model selection criteria (pure ML, AIC, AICc, BIC, HQC) for the 
  different sub-classes, model selection is done for each of the criteria and 
  the parameter estimates for each distribution sub-class are returned.
}
\usage{
pearsonMSC(x, ...)
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
  For the ML estimation, see the details of \code{\link{pearsonFitML}}.
  The considered Model Selection Criteria (MSCs) are 'pure' Maximum Likelihood
  (\code{ML}), Akaike Information Criterion (\code{AIC}), corrected AIC 
  (\code{AICc}), Bayes Information Criterion (\code{BIC}, also known as Schwarz 
  Criterion), and Hannan-Quinn-Criterion (\code{HQC}). The definitions used 
  for the different MSCs are
  \itemize{
    \item for \code{ML}: \eqn{-2\cdot \ln L(\theta)}{-2 ln L(\theta)}
    \item for \code{AIC}: \eqn{-2\cdot \ln L(\theta)+2\cdot k}{-2 ln L(\theta) + 2 k}
    \item for \code{AICc}: \eqn{-2\cdot \ln L(\theta)+2\cdot k\cdot\frac{n}{n-k-1}}{-2 ln L(\theta) + 2 k n / (n-k-1)}
    \item for \code{BIC}: \eqn{-2\cdot \ln L(\theta)+ k\cdot \ln(n)}{-2 ln L(\theta) + k ln(n)}
    \item for \code{HQC}: \eqn{-2\cdot \ln L(\theta)+2\cdot k\cdot \ln(\ln(n))}{-2 ln L(\theta) + 2 k ln(ln(n))}
  }
  where \eqn{\ln L(\theta)}{ln L(\theta)} denotes the log-Likelihood,
  \eqn{n}{n} denotes the number of observations (ie, the length of \code{x}) 
  and \eqn{k}{k} denotes the number of parameters of the distribution 
  (sub-class).
  
  The best model \emph{minimizes} the corresponding MSC function values.
}
\value{
A list containing
  \item{MSCs }{a matrix with rows \code{ML}, \code{AIC}, \code{AICc}, 
    \code{SIC}, \code{HQC} for the different model selection criteria and 
    columns \code{0}, \code{I}, \code{II}, \code{III}, \code{IV}, \code{V}, 
    \code{VI},\code{VII} for the different distribution types containing the 
    values of the criterion/distribution type-combinations.}
  \item{logLik }{a vector with the log-likelihood values for 
    the different distribution types.}
  \item{FittedDistributions}{a list with vectors of the parameter estimates 
    (preceeded by the distribution type number) for the 8 Pearson distribution
    sub-classes.}
  \item{Best}{a list with components \code{ML}, \code{AIC}, \code{AICc},
    \code{SIC}, and \code{HQC} containing the 'best' distributions 
    (distribution type number and parameter estimates) for the different 
    Model Selection Criteria.}
}
%%\references{
%%  TBD
%%}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de}
}
\note{
  The implementation is still preliminary (and slow). No analytical results are 
  used, ie. no analytical solutions for ML estimators and no analytical 
  gradients.
}\seealso{
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}},
  \code{\link[PearsonDS]{pearsonFitML}}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## Generate sample 
DATA <- rpearson(1000,moments=c(mean=1,variance=2,skewness=1,kurtosis=5))
## Call pearsonMSC for model selection
MSC  <- pearsonMSC(DATA,control=list(iter.max=1e5,eval.max=1e5))
## log-Likelihood values for all distribution sub-classes
print(MSC$logLik)
## Values for all MSCs and distribution sub-classes
print(MSC$MSCs)
## Model selection for all MSCs
print(MSC$Best)
}
\keyword{ distribution }
