\name{apc}
\alias{apc}
\title{All Possible Comparisons Analysis}
\description{apc() applies the all possible comparisons procedure to identify the active effects.}
\usage{
apc(y, x, maxsize, level=0.05, method=1, data=NULL, effnames=NULL, reps=50000, dp=4)
}

\arguments{
  \item{y}{Either the response vector or the model formula for the full model.} 
  \item{x}{The model matrix for the full model - only used  when y is a response vector.}
  \item{maxsize}{The maximum model size.}
  \item{level}{The level of error control.}
  \item{method}{The type of error control: 1 = IER, 2 = EER, 3 = FDR.}
  \item{data}{Optional data frame}
  \item{effnames}{Optional vector containing labels for the candidate effects.}
  \item{reps}{The number of repetitions used by the Monte Carlo simulation algorithm which estimates the set of penalties (default is 50000).}
  \item{dp}{the number of decimal places returned for estimates of effects. Default is 4.}
}
\details{
The APC procedure is based on minimizing an AIC-like criterion: APC = log(ResSS) + p where p is a penalty term that increases as the size of the candidate model increases. The penalties can be selected to control either the individual error rate (IER), the experimentwise error rate (EER) or the false discovery rate (FDR) at a specified level. In addition to the type and level of error control, the penalties also depend on the run size of the experiment, the number of candidate effects and the maximum model size.
}
\value{
A list with components
\item{Results}{A data frame that summarizes the results of the APC analysis. The best model of each size 
      is indicated along with its ResSS and value for the APC criterion.}
\item{Penalties}{A vector containing the penalties used for the APC procedure.}
\item{level}{The level of error control.}
\item{ErrorType}{The type of error control used.}
\item{k}{The number of candidate effects.}
\item{m}{The maximum model size.}
\item{apc}{The value of APC for the selected model}
\item{Ests}{A vector containing the estimated effects}
\item{ActEffs}{A vector containing the names of the effects included in the selected model.}
\item{NonActEffs}{A vector containing the names of the effects not included in the selected model.}
}
\references{
Miller A.E. (2005) _The analysis of unreplicated factorial
experiments using all possible comparisons_
Technometrics, 47, 51-63, 2005.
}

\author{Arden Miller and Abu Zar Md. Shafiullah}

\note{
Penalties are estimated using Monte Carlo simulations and thus the estimates will not be exactly the same each time the function is run. The precision of the estimates can be increased by increasing the number of reps but the function will take longer to run.  The amount of time needed to run this programme increases as the values of n, k and m increase. For larger experiments it may be necessary to reduce the number of reps.
}


\seealso{
\code{\link{IERpenalties}}, \code{\link{EERpenalties}}, \code{\link{FDRpenalties}}.
}

\examples{
## This example demonstrates the analysis of an artificial data set for a unreplicated 
## factorial design for four two-level factors. The values of "resp" were generated as
## "resp<-round(10+8*x1+5*x3+7*x4+6*x1*x4+rnorm(16), 2)". The data is contained in the 
## data frame "testdata". A maximum model size of 6 and an IER of .05 are used.
apc(resp~x1*x2*x3*x4, maxsize=6, data=testdata, method=1, level=.05, reps=9000)
}

\keyword{design}
\keyword{models}
