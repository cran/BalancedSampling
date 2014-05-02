\name{lcube}
\alias{lcube}
\docType{package}
\title{
Local cube method (Doubly balanced sampling)
}
\description{
Select doubly balanced samples with prescribed inclusion probabilities from a finite population. To have a fixed sample size, include the inclusion probabilities as a balancing variable in \code{Xbal} and make sure the inclusion probabilities sum to a positive integer. This is a simplified (optimized for speed) implementation of the local cube method (doubly balanced sampling). Landing is done by dropping balancing variables (from rightmost column, so keep inclusion probabilities in first column to guarantee fixed size). Euclidean distance is used in the \code{Xspread} space.
}
\usage{
lcube(prob,Xspread,Xbal)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{Xspread}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and r columns}
}
\value{
Returns a vector of selected indexes in 1,2,...,N. 
}
\references{
Grafström, A. and Tillé, Y. (2013). Doubly balanced spatial sampling with spreading and restitution of auxiliary totals. Environmetrics, 24(2), 120-131.
}
\examples{
\dontrun{
# Example 1
set.seed(12345);
N = 1000; # population size
n = 100; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(runif(N),runif(N)); # matrix of auxiliary variables
s = lcube(p,X,cbind(p)); # select sample 
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
  s = lcube(p,X,cbind(p));
  ep[s]=ep[s]+1;
}
print(ep/nrs);

}
}
