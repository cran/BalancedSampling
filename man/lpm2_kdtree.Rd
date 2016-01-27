\name{lpm2_kdtree}
\alias{lpm2_kdtree}
\title{
Local Pivotal Method 
}
\description{
  The local pivotal method provides a way to perform balanced sampling.  This
  implementation replace linear searches in lpm2, with k-d trees.  K-d trees
  are binary trees used to effectively search high dimensional spaces, and 
  reduce the average computational complexity of lpm2 from O(N^2) to 
  O(N log(N)).
}
\usage{

  lpm2_kdtree(
     prob,
     x, 
     m=40
     ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{
    An array of length N such that the sum of prob is equal to the sample size,
    where the N is the number of rows of x.
  }
  \item{x}{
    A matrix of N rows and p columns, each row is assumed to be a sampling unit. 
  }
  \item{m}{
    Max leaf size used as a terminal condition for building the k-d tree.
  }
}
\value{ 
Returns a vector of selected indexes from the matrix x.  The results are 
identical to the lpm2 function when no ties in distance occur.
}
\references{
Lisic, J., Jonathan (2015)  Parcel Level Agricultural Land Cover Prediction. 
(Unpublished doctoral dissertation). George Mason University, Fairfax, Virginia.  
}
\author{
  Jonathan Lisic
}
\examples{


N <- 1000
n <- 100
x <- cbind( runif(N), runif(N)) 

set.seed(100)
Cprog <- proc.time()
sampled <- lpm2_kdtree( rep(n/N,N), x)
print("lpm2_kdtree running time")
print(proc.time() - Cprog) 

}
