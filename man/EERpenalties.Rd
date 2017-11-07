\name{EERpenalties}
\alias{EERpenalties}
\title{EER Penalties for APC}

\description{
EERpenalties() generates the set of penalties used for the APC criterion so that the experimentwise error rate (EER) is controlled at the specified level. The values of the penalties are estimated using a Monte Carlo simulation procedure and also depend on the run size of the experiment, the number of candidate effects and the maximum model size.
} 

\usage{
EERpenalties(n, k = n - 1, m = min(n - 2, k), eer = 0.2, reps = 50000, rnd = 3)
}

\arguments{
\item{n}{The number of experimental runs in the study, i.e. the row dimension of the design matrix (for orthogonal 2-level designs n must be a multiple of 4).}
\item{k}{The number of candidate effects under study, i.e., the column dimension of the design matrix (1 < k <= n-1).}
\item{m}{The maximum size of the candidate models (0 < m < min(n-2, k)).}
\item{eer}{The level (0 < eer < 1) at which the experimentwise error rate will be controlled (default is eer = 0.2).}
\item{reps}{The number of repetitions used by the Monte Carlo simulation algorithm which estimates the set of penalties (default is 50000).}
\item{rnd}{The number of decimal places returned for the estimated penalties (default is rnd = 3).}
}

\value{
A vector of containing the m + 1 penalties for the APC procedure that controls the EER at the specified level.
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
\code{\link{apc}}, \code{\link{IERpenalties}}, \code{\link{FDRpenalties}}
}
\examples{
## Penalties for a 8-run experiment that has 5 candidate effects are generated.
## The maximum model size is set to 5 and an experimentwise error rate of .25 is used.
EERpenalties(n = 8, k = 5, m = 5, eer = .25, reps = 15000)
}

\keyword{ design }
\keyword{ models }
