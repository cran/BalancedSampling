\name{lcps}
\alias{lcps}
\docType{package}
\title{Locally correlated Poisson samling}
\description{
  Selects spatially balanced sampling with prescribed inclusion
  probabilites from a finite population using Locally correlated
  Poisson sampling, a variant of SCPS (\code{\link[=scps]{scps()}}).
}
\usage{
  lcps(prob, x)
}
\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{x}{matrix of (standardized) auxiliary variables of N rows and q columns}
}
\value{
  Returns a vector of selected indices in 1,2,...,N.
  If the inclusion probabilities sum to n, where n is integer, then the
  sample size is fixed (n).
}
\details{
\code{lcps} uses euclidean distance on auxiiary variables to calculate distance between units.
}
\references{
  Prentius, W. (2022). Locally correlated Poisson sampling. Manuscript.
}
\examples{
\dontrun{
  N = 1000L; # size of the population
  prob = runif(N, 0.1, 0.3); # Inclusion probabilities
  x = matrix(rnorm(3*N), ncol = 3L); # N * 3 matrix of auxiliary variables
  s = lcps(prob, x); # sample indices
}
}
