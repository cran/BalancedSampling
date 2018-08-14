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
  O(N log(N)).  Both nearest neighbor and approximate nearest neighbor 
  searching algorithms are provided.
}
\usage{

  lpm2_kdtree(
     prob,
     x, 
     m=40,
     algorithm = "kdtree", 
     maxCheck = 4,
     termDist = 0.1
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
  \item{algorithm}{
    The algorithm used to search the k-d tree.  The algorithms include "kdtree",
    "kdtree-count", and "kdtree-dist".  The "kdtree" algorithm reproduces the
   lpm2 using a k-d tree for nearest neighbor search.  "kdtree-count" and 
   "kdtree-dist" use approximate nearest neighbor searches based on number of
   nodes to check and minimal sufficient distance respectfully.  
  }
  \item{maxCheck}{
    A positive integer scalar parameter only used when the algorithm "kdtree-count"
    is specified.  This parameter is the maximum number of non-empty leaf nodes 
    to check for a nearest neighbor.
  }
  \item{termDist}{
    A positive valued scalar parameter only used when the algorithm "kdtree-dist"
    is specified.  This parameter specifies a minimal sufficient distance to be 
    considered a nearest neighbor.  No tie handling is performed; the first 
    nearest neighbor found will be used. 
  }
}
\value{ 
Returns a vector of selected indexes from the matrix x.  If the algorithm is
"kdtree", the results are identical to the lpm2 function when no ties in distance 
occur.
}
\references{
Lisic, J., Jonathan (2015)  Parcel Level Agricultural Land Cover Prediction. 
(Unpublished doctoral dissertation). George Mason University, Fairfax, Virginia.  
}
\author{
  Jonathan Lisic
}
\examples{
\dontrun{

N <- 1000
n <- 100
x <- cbind( runif(N), runif(N)) 

set.seed(100)
Cprog <- proc.time()
sampled <- lpm2_kdtree( rep(n/N,N), x)
print("lpm2_kdtree running time")
print(proc.time() - Cprog) 

}
}
