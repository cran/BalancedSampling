\name{cubestratified}
\alias{cubestratified}
\docType{package}
\title{
Stratified balanced sampling with pooling of landing phases
}
\description{
This is a fast implementation of stratified balanced sampling. To have a fixed sample size, include the inclusion probabilities as a balancing variable in \code{Xbal} and make sure the inclusion probabilities sum to a positive integer (within each stratum). 
}
\usage{
cubestratified(prob,Xbal,integerStrata)	
}

\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{Xbal}{matrix of balancing auxiliary variables of N rows and r columns}
  \item{integerStrata}{vector of length N with stratum number}
}
\value{
Returns a vector of length N with sampling indicators. 
}
\references{
Chauvet, G. (2009). Stratified balanced sampling. Survey Methodology, 35, 115-119.
}
\examples{
\dontrun{
# Example 1
N = 10;
n = 5;
p = rep(n/N,N);
strata = c(1,1,2,2,3,3,4,4,5,5);
indicators = cubestratified(p,cbind(p),strata);
s = (1:N)[indicators==1];
}
}
