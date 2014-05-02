\name{lcubeflightphase}
\alias{lcubeflightphase}
\docType{package}
\title{
Flight phase for the local cube method 
}
\description{
Flight phase for the local cube method. To have a fixed sample size, include the inclusion probabilities as a balancing variable in \code{Xbal} and make sure the inclusion probabilities sum to a positive integer. This is a simplified (optimized for speed) implementation of the flight phase of the local cube method (doubly balanced sampling). Euclidean distance is used in the \code{Xspread} space.
}
\usage{
lcubeflightphase(prob,Xspread,Xbal)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{Xspread}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and r columns}
}
\value{
Returns a vector of length N with new probabilities, where at most r are non-integer. 
}
\references{
Grafström, A. and Tillé, Y. (2013). Doubly balanced spatial sampling with spreading and restitution of auxiliary totals. Environmetrics, 24(2), 120-131.
}
\examples{
\dontrun{
# Example 1
# Select sample
set.seed(12345);
N = 100; # population size
n = 10; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(runif(N),runif(N)); # matrix of auxiliary variables
pflight = lcubeflightphase(p,X,cbind(p,X)); 
# check balance
colSums(X)
colSums(X/p*pflight)

}
}
