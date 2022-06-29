\name{scps_getrand}
\alias{scps_getrand}
\docType{package}
\title{
Spatially correlated Poisson sampling
}
\description{
Retrieves a vector of random numbers resulting in the provided scps sample
}
\usage{
scps_getrand(prob,x,s)	
}
\arguments{
  \item{prob}{vector of length N with inclusion probabilities}
  \item{x}{matrix of (standardized) auxiliary variables of N rows and q columns}
  \item{s}{vector of length N with inclusion indicators}
}
\value{
Retrieves a vector of length N with random numbers resulting in the provided scps sample
}
\references{
Grafström, A. (2012). Spatially correlated Poisson sampling. Journal of Statistical Planning and Inference, 142(1), 139-147.
}
