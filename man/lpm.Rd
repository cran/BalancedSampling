\name{lpm}
\alias{lpm}
\docType{package}
\title{
Local pivotal method (sub-optimal)
}
\description{
Select spatially balanced samples with prescribed inclusion probabilities from a finite (large) population using a sub-optimal implementation of the local pivotal method. Euclidean distance is used in the \code{x} space.
}
\usage{
lpm(prob,x,h)	
}
\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{x}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{h}{positive integer, size of window in the list to search for nearest neighbor}
}
\value{
Returns a vector of selected indexes in 1,2,...,N. If the inclusion probabilities sum to n, where n is integer, then the sample size is fixed (n).
}
\examples{
\dontrun{
# Example 1
set.seed(12345);
N = 1000; # population size
n = 100; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(runif(N),runif(N)); # matrix of auxiliary variables
h = 100; # size of search window (for finding nearest neighbor)
s = lpm(p,X,h); # select sample 
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
  s = lpm(p,X,10);
  ep[s]=ep[s]+1;
}
print(ep/nrs);
}
}
