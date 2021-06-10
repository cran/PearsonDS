\name{PearsonIV}
\Rdversion{1.1}
\alias{PearsonIV}
\alias{dpearsonIV}
\alias{ppearsonIV}
\alias{qpearsonIV}
\alias{rpearsonIV}
\title{
  The Pearson Type IV Distribution
}
\description{
  Density, distribution function, quantile function and random generation for 
  the Pearson type IV distribution.
}
\usage{

dpearsonIV(x, m, nu, location, scale, params, log = FALSE)

ppearsonIV(q, m, nu, location, scale, params, lower.tail = TRUE, 
           log.p = FALSE, tol = 1e-08, ...)

qpearsonIV(p, m, nu, location, scale, params, lower.tail = TRUE, 
           log.p = FALSE, tol = 1e-08, ...)

rpearsonIV(n, m, nu, location, scale, params)
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
  \item{m}{
    first shape parameter of Pearson type IV distribution.
}
  \item{nu}{
    second shape parameter (skewness) of Pearson type IV distribution.
}
  \item{location}{
    location parameter of Pearson type IV distribution.
}
  \item{scale}{
    scale parameter of Pearson type IV distribution.
}
  \item{params}{
    vector/list of length 4 containing parameters \code{m}, \code{nu}, 
    \code{location}, \code{scale} for Pearson type IV distribution 
    (in this order!). 
}
  \item{log, log.p}{
    logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
}
  \item{lower.tail}{
    logical; if \code{TRUE}, probabilities are \eqn{P[X\le x]}{P[X<=x]}, 
    otherwise, \eqn{P[X>x]}{P[X>x]}.
}
  \item{tol}{
    relative tolerance for evaluation of hypergeometric function 2F1 
    (\code{ppearsonIV}) or absolute target \code{q}-error for Newton method 
    (\code{qpearsonIV}).
}
  \item{\dots}{
    further parameters for underlying hypergeometric function.
}
}
\details{
  The Pearson type IV distribution with location parameter 
  \code{location}\eqn{=\lambda}{=lamda},
  scale parameter \code{scale}\eqn{=a}, and shape parameters 
  \eqn{m} and \eqn{\nu}{nu} can
  be obtained by its probability density function
  \deqn{f(x) = \frac{\left|\frac{\Gamma(m+\frac{\nu}{2}i)}{\Gamma(m)}\right|^2}
        {a B(m-\frac{1}{2},\frac{1}{2})}
        \left[1+\left(\frac{x-\lambda}{a}\right)^2\right]^{-m} 
        e^{-\nu \mathrm{arctan}\left(\frac{x-\lambda}{a} \right)}}{f(x) =
        |Gamma(m+nu/2*i)/Gamma(m)|^2/(a*B(m-1/2,1/2))*[1+((x-lambda)/a)^2]^(-m) 
        * exp[-nu*arctan((x-lambda)/a)]}
  for \eqn{a>0}, \eqn{m>\frac{1}{2}}{m>1/2}, \eqn{\nu\ne 0}{nu<>0} 
  (\eqn{\nu=0}{nu=0} corresponds to the Pearson type VII distribution family).
  
  The normalizing constant, which involves the complex Gamma function, is
  calculated with \code{lngamma_complex} (see \code{\link[gsl]{Gamma}}) of
  package \code{gsl}, if package \code{gsl} is installed.  
  Section 5.1 of [2] contains an algorithm (C code) for the calculation of the
  normalizing constant, which is used otherwise, but this will be very slow for 
  large absolute values of \eqn{\nu}{nu}.
  
  The generation of random numbers (\code{rpearsonIV}) uses the C code from 
  section 7 of [2]. It is (thus) restricted to distributions with \eqn{m>1}.
  
  For the cumulative distribution function (\code{ppearsonIV}), numerical 
  integration of the density function is used, if package \code{gsl} is not
  available. 
  If package \code{gsl} is installed, three different
  methods are used, depending on the parameter constellation (the corresponding 
  parameter regions were obtained by comprehensive benchmarks):
  \itemize{
    \item numerical integration of the density function
    \item cdf representation of Heinrich [2]
    \item cdf representation of Willink [4]
  }
  The hypergeometric functions involved in the latter two representations are 
  approximated 
  via partial sums of the corresponding series (see [1], 15.1.1, p. 556). 
  Depending on the parameter constellation, transformation 15.3.5 of [1] 
  (p. 559) is applied for Heinrich's method.
  The evaluation of the partial sums is first carried out in (ordinary) 
  double arithmetic. If cancellation reduces accuracy beyond \code{tol}, the
  evaluation is redone in double-double arithmetics. If cancellation still
  reduces accuracy beyond \code{tol}, the evaluation is again redone in 
  quad-double arithmetic. Code for double-double and quad-double arithmetics 
  is based on [3].
  For Willink's representation, the hypergeometric function in the denominator 
  of \eqn{R} in equation (10) is evaluated via complex gamma 
  functions (see [1], 15.1.20, p. 556), which is fast and much more stable. 
  A warning is issued if the approximation of the hypergeometric function 
  seems to fail (which should not 
  happen, since numerical integration should be carried out for critical
  parameter constellations).
  
  The quantile function (\code{qpearsonIV}) is obtained via Newton's method.
  The Newton iteration begins in the (single) mode of the distribution, 
  which is easily calculated (see [2], section 3). Since the mode of the 
  distribution is the only inflection point of the cumulative distribution
  function, convergence is guaranteed.
  The Newton iteration stops when the target \code{q}-error is reached 
  (or after a maximum of 30 iterations).
}
\value{
  \code{dpearsonIV} gives the density, \code{ppearsonIV} gives the distribution 
  function, \code{qpearsonIV} gives the quantile function, and \code{rpearsonIV} 
  generates random deviates.
}
\note{
  Many calculations are done in logarithms to avoid IEEE 754 underflow/overflow
  issues.
  
  The description of quad-double arithmetics in [3] contains minor errors:
  in algorithm 9 (p. 6), lines 9 and 10 should be interchanged; 
  in algorithm 14 (p. 10), \eqn{k<2} should be replaced with \eqn{k<3} (line 10)
  and \eqn{k<3} should be replaced with \eqn{k<4} (line 11).
}
\section{Warning}{  
  If at all possible, package \code{gsl} should be installed. 
  Otherwise, calculations for distributions with (more or less) extreme 
  parameters may slow down by factors of more than 1000 and/or may become
  unstable.
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

  [4] Willink, R. (2008) \emph{A Closed-form Expression for the Pearson Type IV 
  Distribution Function}, Aust. N. Z. J. Stat. 50 (2), pp. 199-205
}
\author{
Martin Becker \email{martin.becker@mx.uni-saarland.de},
Stefan Klößner \email{S.Kloessner@mx.uni-saarland.de} and
Joel Heinrich
}
\seealso{
  \code{\link{PearsonDS-package}},
  \code{\link[PearsonDS]{Pearson}}
}
\examples{
## define Pearson type IV parameter set with m=5.1, nu=3, location=0.5, scale=2
pIVpars <- list(m=5.1, nu=3, location=0.5, scale=2)
## calculate probability density function
dpearsonIV(-2:2,params=pIVpars)
## calculate cumulative distribution function
ppearsonIV(-2:2,params=pIVpars)
## calculate quantile function
qpearsonIV(seq(0.1,0.9,by=0.2),params=pIVpars)
## generate random numbers
rpearsonIV(5,params=pIVpars)
}
\keyword{ distribution }
