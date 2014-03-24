\name{sb}
\alias{sb}
\docType{package}
\title{
Spatial balance
}
\description{
Calculates spatial balance of a sample subject to inclusion probabilities and auxiliary space
}
\usage{
sb(p,x,s)	
}
\arguments{
  \item{p}{vector of length N with inclusion probabilities}
  \item{x}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{s}{the sample, vector of length n}
}
\value{
Number, the spatial balance
}
\references{
Grafström, A., Lundström, N.L.P. and Schelin, L. (2012). Spatially balanced sampling through the Pivotal method. Biometrics 68(2), 514-520.
}
\examples{
\dontrun{
# check spatial balance
set.seed(1234567);
N = 500;
n = 70;
p = rep(n/N,N);
X = cbind(runif(N),runif(N));

# select lpm1 sample
s = lpm1(p,X);
# calculate balance
B = sb(p,X,s);

}
}
