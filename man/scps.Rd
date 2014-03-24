\name{scps}
\alias{scps}
\docType{package}
\title{
Spatially correlated Poisson sampling
}
\description{
Select spatially balanced samples with prescribed inclusion probabilities from a finite population. This implementation uses the maximal weight strategy and Euclidean distance.
}
\usage{
scps(prob,x)	
}
\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{x}{matrix of (standardized) auxiliary variables of N rows and q columns}
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
s = scps(p,X); # select sample 
plot(X[,1],X[,2]); # plot population
points(X[s,1],X[s,2], pch=19); # plot sample

# Example 2
# check inclusion probabilities
set.seed(12345);
p = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9); # prescribed inclusion probabilities
N = length(p); # population size
X = cbind(runif(N),runif(N)); # some artificial auxiliary variables
ep = rep(0,N); # empirical inclusion probabilities
nrs = 10000; # repetitions
for(i in 1:nrs){
  s = scps(p,X);
  ep[s]=ep[s]+1;
}
print(ep/nrs);
}
}
