\name{lcubelandingphase}
\alias{lcubelandingphase}
\docType{package}
\title{
Landing phase for the local cube method 
}
\description{
Landing is done by dropping balancing variables (from rightmost column). Euclidean distance is used in the \code{Xspread} space.
}
\usage{
lcubelandingphase(prob,probflight,Xspread,Xbal)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{probflight}{vector of length N with probabilities from flightphase}
  \item{Xspread}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and r columns}
}
\value{
Returns a vector of length N with indicators. 
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
pflight = lcubeflightphase(p,X,cbind(p,X)); # flight
indicators = lcubelandingphase(p,pflight,X,cbind(p,X)); # landing
# final sample
s = (1:N)[indicators==1];
}
}
