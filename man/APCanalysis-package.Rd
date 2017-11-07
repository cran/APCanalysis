\docType{package}
\name{APCanalysis-package}
\alias{APCanalysis-package}
\alias{APCanalysis}
\title{Analysis of Unreplicated Orthogonal Experiments using All Possible Comparisons} 
\description{
This package provides functions to analyse data from unreplicated orthogonal experiments such as 2-level factorial and fractional factorial designs and Plackett-Burman designs using the all possible comparisons (APC) methodology.
}
\details{
apc() identifies the active effects from an unreplicated orthogonal experiment using a modified version of the all possible comparisons (APC) procedure proposed by Miller (2005). This function has been designed specifically to analyse data from two-level designs including full factorial designs, regular fractional factorial designs and Plackett-Burman designs. 

The APC procedure is based on minimizing an AIC-like criterion: APC = log(ResSS) + p where p is a penalty term that increases as the size of the candidate model increases. The APC procedure can be adapted to control either the individual error rate (IER), the experimentwise error rate (EER) or the false discovery rate (FDR). The functions IERpenalties(), EERpenalties() and FDRpenalties() can be used to estimate the penalties used in the APC criterion for each type of error control.
}
\author{
Arden Miller and Abu Zar Md. Shafiullah.

Maintainer: Arden Miller <a.miller@auckland.ac.nz>
}
\references{
Miller A.E. (2005) _The analysis of unreplicated factorial
experiments using all possible comparisons_
Technometrics, 47, 51-63, 2005.
}

\keyword{ package }
\keyword{ design }
\seealso{
\code{\link{apc}}, \code{\link{IERpenalties}}, \code{\link{EERpenalties}}, \code{\link{FDRpenalties}}.
}
\examples{
## This example demonstrates the analysis of an artificial data set for a  12-run 
## Plackett-Burman design stored in "PB12matrix". The values of "PB12response" were 
## generated using the following active effects: B=7, D=5, H=11, I=4 and K=6. The 
## remaining columns were all set to be inactive (effects equal 0). 
my.apc = apc(PB12response, PB12matrix, maxsize=6, method = 2, level = 0.20, reps = 10000)
summary(my.apc)
plot(my.apc)
}
