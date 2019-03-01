\name{scps_coord}
\alias{scps_coord}
\docType{package}
\title{
Spatially correlated Poisson sampling with coordination
}
\description{
Select spatially balanced samples with prescribed inclusion probabilities from a finite population. This implementation uses the maximal weight strategy and Euclidean distance. To be used with permanent random numbers.
}
\usage{
scps_coord(prob,x,rand)	
}
\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{x}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{rand}{vector of length N with permanent random numbers}
}
\value{
Returns a vector of selected indexes in 1,2,...,N. If the inclusion probabilities sum to n, where n is integer, then the sample size is fixed (n).
}
\references{
Grafström, A. (2012). Spatially correlated Poisson sampling. Journal of Statistical Planning and Inference, 142(1), 139-147.
}
\examples{
\dontrun{
# Example 1
set.seed(12345);
N = 1000; # population size
n = 100; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(runif(N),runif(N)); # matrix of auxiliary variables
u = runif(N);
s = scps_coord(p,X,u); # select sample 
plot(X[,1],X[,2]); # plot population
points(X[s,1],X[s,2], pch=19); # plot sample
}
}
