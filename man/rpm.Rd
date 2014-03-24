\name{rpm}
\alias{rpm}
\docType{package}
\title{
Random pivotal method
}
\description{
Select samples with prescribed inclusion probabilities from a finite population. This design has high entropy. In each of the (at most) N steps, two undecided units are selected at random to compete.
}
\usage{
rpm(prob)	
}
\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
}
\value{
Returns a vector of selected indexes in 1,2,...,N. If the inclusion probabilities sum to n, where n is integer, then the sample size is fixed (n).
}
\examples{
\dontrun{
# Example 1
set.seed(12345);
N = 100; # population size
n = 10; # sample size
p = rep(n/N,N); # inclusion probabilities
s = rpm(p); # select sample 

# Example 2
# check inclusion probabilities
set.seed(12345);
p = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9); # prescribed inclusion probabilities
N = length(p); # population size
ep = rep(0,N); # empirical inclusion probabilities
nrs = 10000; # repetitions
for(i in 1:nrs){
  s = rpm(p);
  ep[s]=ep[s]+1;
}
print(ep/nrs);
}
}
