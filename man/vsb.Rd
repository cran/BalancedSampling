\name{vsb}
\alias{vsb}
\docType{package}
\title{
Variance estimator for spatially balanced sample
}
\description{
Variance estimator of HT estimator of population total of target variable y
}
\usage{
vsb(probs,ys,xs)	
}
\arguments{
  \item{probs}{vector of length n (sample) with inclusion probabilities}
  \item{ys}{vector of target variable y of length n (sample)}
  \item{xs}{matrix of (standardized) auxiliary variables of n rows (sample) and q columns}
}
\value{
Number, the estimated variance
}
\references{
Grafström, A., and Schelin, L. (2014). How to select representative samples. Scandinavian Journal of Statistics.
}
