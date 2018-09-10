\name{hlpm}
\alias{hlpm}
\docType{package}
\title{
Hierarchical local pivotal method
}
\description{
Hierarchical local pivotal method (hlpm) selects an initial sample using the local pivotal method and then splits the sample into subsamples of given sizes using a successive (hierarchical) selection with the local pivotal method. Can be used with any prescribed inclusion probabilities that sum to an integer n. The sizes of the subsamples must also sum to n. It is used to select several subsamples such that each subsample is spatially balanced and the combined sample (the union of the subsamples) is also spatially balanced. Licence (GPL >=2).
}
\usage{
hlpm(p,X,sizes)
}
\arguments{
  \item{p}{vector of inclusion probabilities for initial sample.}
  \item{X}{matrix of auxiliary variables.}
  \item{sizes}{vector of sizes of subsamples whose sum must match the sum of the initial inclusion probabilities.}
}
\value{Returns a list with population indexes of initial sample S and a vector sampleNumber indicating the number of the subsample of each unit.}

\examples{
\dontrun{
############
## Example with two subsamples
############

N = 100; # population size
X = cbind(runif(N),runif(N)); # auxiliary variables
n = 10; # size of initial sample
p = rep(n/N,N); # inclusion probabilities of initial sample
sizes = c(7,3); # sizes of the two subsamples
hlpm(p,X,sizes) # selection of samples using hierarchical local pivotal method

}
}
