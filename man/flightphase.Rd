\name{flightphase}
\alias{flightphase}
\docType{package}
\title{
Flight phase of the cube method
}
\description{
This is a fast implementation of the flight phase of the cube method. To have a fixed sample size, include the inclusion probabilities as a balancing variable in \code{Xbal} and make sure the inclusion probabilities sum to a positive integer. 
}
\usage{
flightphase(prob,Xbal)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and q columns}
}
\value{
Returns a vector of length N with new probabilities, where at most q are non-integer. 
}
\references{
Deville, J. C. and  Tillé, Y. (2004). Efficient balanced sampling: the cube method. Biometrika, 91(4), 893-912.

Chauvet, G. and Tillé, Y. (2006). A fast algorithm for balanced sampling. Computational Statistics, 21(1), 53-62.
}
\examples{
\dontrun{
# Example 1
# Select sample and check balance
set.seed(12345);
N = 100; # population size
n = 10; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(p,runif(N),runif(N)); # matrix of auxiliary variables

pflight = flightphase(p,X);

# check balance after flight
colSums(X)
colSums(X/p*pflight)

# select final sample as indicators
indicators = landingphase(p,pflight,X);

# check final balance
colSums(X)
colSums(X/p*indicators)

# final sample as indexes
s = (1:N)[indicators==1];
}
}
