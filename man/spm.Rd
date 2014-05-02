\name{spm}
\alias{spm}
\docType{package}
\title{
Sequential pivotal method (also known as ordered pivotal sampling and Deville's systematic sampling)
}
\description{
Select samples with prescribed inclusion probabilities from a finite population. The resulting samples are well spread in the list (similar to systematic sampling). In each of the (at most) N steps, two undecided units with smallest index are selected to compete.
}
\usage{
spm(prob)	
}
\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
}
\value{
Returns a vector of selected indexes in 1,2,...,N. If the inclusion probabilities sum to n, where n is integer, then the sample size is fixed (n).
}
\references{
Deville, J.-C. and  Tillé, Y. (1998). Unequal probability sampling without replacement through a splitting method. Biometrika 85, 89-101.

Chauvet, G. (2012). On a characterization of ordered pivotal sampling. Bernoulli, 18(4), 1320-1340.
}
\examples{
\dontrun{
# Example 1
set.seed(12345);
N = 100; # population size
n = 10; # sample size
p = rep(n/N,N); # inclusion probabilities
s = spm(p); # select sample 

# Example 2
# check inclusion probabilities
set.seed(12345);
p = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9); # prescribed inclusion probabilities
N = length(p); # population size
ep = rep(0,N); # empirical inclusion probabilities
nrs = 10000; # repetitions
for(i in 1:nrs){
  s = spm(p);
  ep[s]=ep[s]+1;
}
print(ep/nrs);
}
}
