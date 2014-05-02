\name{landingphase}
\alias{landingphase}
\docType{package}
\title{
Landing phase of the cube method
}
\description{
Landing is done by dropping balancing variables (from rightmost column).  
}
\usage{
landingphase(prob,probflight,Xbal)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{probflight}{vector of length N obtained from the flightphase}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and q columns}
}
\value{
Returns a vector of length N with inclusion indicators. 
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
N = 100; # population size
n = 10; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(p,runif(N),runif(N)); # matrix of auxiliary variables
pflight = flightphase(p,X); # flight
indicators = landingphase(p,pflight,X); # landing
# final sample
s = (1:N)[indicators==1];
}
}
