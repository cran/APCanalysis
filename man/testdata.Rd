\name{testdata}
\alias{testdata}
\docType{data}
\title{ testdata data frame }
\description{
A constructed data frame to illustrate the use of the functions in the APC package.
}
\usage{data("testdata")}
\format{
The "testdata" data frame has 16 rows and 5 columns:
\itemize{
\item{resp} response variable.
\item{x1} binary (-1 or +1) explanatory variable 1.
\item{x2} binary (-1 or +1) explanatory variable 2.
\item{x3} binary (-1 or +1) explanatory variable 3.
\item{x4} binary (-1 or +1) explanatory variable 4.
}
}
\details{
The values of "resp" were generated as
resp = round(10+8*x1+5*x3+7*x4+6*x1*x4+rnorm(16), 2).`
}
\examples{
data(testdata)
}
\keyword{datasets}

