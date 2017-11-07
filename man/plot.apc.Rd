\name{plot.apc}
\alias{plot.apc}
\title{Graphical Summary of an APC Analysis}
\description{Produces a scatterplot of minimum APC versus model size. This is useful for visualizing the relative 
             values of APC for the best models of each size.}
\usage{
\method{plot}{apc}(x, elabs = TRUE, ...)
}

\arguments{
  \item{x}{apc object}
  \item{elabs}{use effect labels as plotting characters}
  \item{\dots}{other arguments}
}

\value{none}
\author{Arden Miller and Abu Zar Md. Shafiullah}

\seealso{
\code{\link{apc}}, \code{\link{summary.apc}}
}
\examples{
## This example demonstrates the analysis of an artificial data set for a unreplicated 
## factorial design for four two-level factors. The values of "resp" were generated as
## "resp<-round(10+8*x1+5*x3+7*x4+6*x1*x4+rnorm(16),2)". The data is contained in the 
## data frame "testdata". A maximum model size of 6 and an IER of .05 are used.
my.apc = apc(resp~x1*x2*x3*x4, maxsize=6, data=testdata, method=1, level=.05, reps=9000)
plot(my.apc)
}

\keyword{design}
\keyword{models}


