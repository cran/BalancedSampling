\name{lcubestratified}
\alias{lcubestratified}
\docType{package}
\title{
Stratified doubly balanced sampling with pooling of landing phases
}
\description{
This is a fast implementation of stratified doubly balanced sampling. To have a fixed sample size, include the inclusion probabilities as a balancing variable in \code{Xbal} and make sure the inclusion probabilities sum to a positive integer (within each stratum). Euclidean distance is used in the \code{Xspread} space.
}
\usage{
lcubestratified(prob,Xspread,Xbal,integerStrata)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{Xspread}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and r columns}
  \item{integerStrata}{vector of length N with stratum number}
}
\value{
Returns a vector of length N with sampling indicators. 
}
\references{
Chauvet, G. (2009). Stratified balanced sampling. Survey Methodology, 35, 115-119.

Grafström, A. and Tillé, Y. (2013). Doubly balanced spatial sampling with spreading and restitution of auxiliary totals. Environmetrics, 24(2), 120-131.
}
\examples{
\dontrun{
# Example 1
N = 10;
n = 5;
p = rep(n/N,N);
Xspread = cbind(1:N);
strata = c(1,1,1,1,1,1,2,2,2,2);
indicators = lcubestratified(p,Xspread,cbind(p),strata);
s = (1:N)[indicators==1];
}
}
