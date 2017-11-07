\name{summary.apc}
\alias{summary.apc}
\title{Summary for an apc object}
\description{Produces a useful summary of an apc object} 
\usage{
\method{summary}{apc}(object, ...)
}

\arguments{
  \item{object}{apc object}
  \item{\dots}{other arguments}
}

\value{none}
\author{Arden Miller and Abu Zar Md. Shafiullah}

\seealso{
\code{\link{apc}}, \code{\link{plot.apc}}
}
\examples{
## This example demonstrates the analysis of an artificial data set for a unreplicated 
## factorial design for four two-level factors. The values of "resp" were generated as
## "resp<-round(10+8*x1+5*x3+7*x4+6*x1*x4+rnorm(16), 2)". The data is contained in the 
## data frame "testdata". A maximum model size of 6 and an IER of .05 are used.
my.apc = apc(resp~x1*x2*x3*x4, maxsize=6, data=testdata, method=1, level=.05, reps=9000)
summary(my.apc)
}

\keyword{design}
\keyword{models}

