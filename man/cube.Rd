\name{cube}
\alias{cube}
\docType{package}
\title{
Cube method (Balanced sampling)
}
\description{
This is a fast implementation of the cube method. To have a fixed sample size, include the inclusion probabilities as a balancing variable in \code{Xbal} and make sure the inclusion probabilities sum to a positive integer. Landing is done by dropping balancing variables (from rightmost column, so keep inclusion probabilities in first column to guarantee fixed size). 
}
\usage{
cube(prob,Xbal)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and r columns}
}
\value{
Returns a vector of selected indexes in 1,2,...,N. 
}
\references{
Deville, J. C. and  Tillé, Y. (2004). Efficient balanced sampling: the cube method. Biometrika, 91(4), 893-912.

Chauvet, G. and Tillé, Y. (2006). A fast algorithm for balanced sampling. Computational Statistics, 21(1), 53-62.
}
\examples{
\dontrun{
# Example 1
# Select sample
set.seed(12345);
N = 1000; # population size
n = 100; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(p,runif(N),runif(N)); # matrix of auxiliary variables
s = cube(p,X); # select sample 


# Example 2
# Check inclusion probabilities
set.seed(12345);
p = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9); # prescribed inclusion probabilities
N = length(p); # population size
ep = rep(0,N); # empirical inclusion probabilities
nrs = 10000; # repetitions
for(i in 1:nrs){
  s = cube(p,cbind(p));
  ep[s]=ep[s]+1;
}
print(ep/nrs);

# Example 3
# How fast is it?
# Let's check with N = 100 000 and 5 balancing variables
set.seed(12345);
N = 100000; # population size
n = 100; # sample size
p = rep(n/N,N); # inclusion probabilities
# matrix of 5 auxiliary variables
X = cbind(p,runif(N),runif(N),runif(N),runif(N)); 
system.time(cube(p,X));  
}
}
