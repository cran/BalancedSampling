\name{inclusionprobabilities}
\alias{inclusionprobabilities}
\docType{package}
\title{
Inclusion probabilities
}
\description{
Computes the first-order inclusion probabilities from a vector of positive numbers (for a probability proportional-to-size sampling design). This function is borrowed from the package "sampling" by Alina Matei and Yves Tillé. Licence (GPL >=2).
}
\usage{
inclusionprobabilities(a,n)
}
\arguments{
  \item{a}{vector of positive numbers}
  \item{n}{sample size}
}

\examples{
\dontrun{
############
## Example 
############
# a vector of positive numbers
a=1:20
# computation of the inclusion probabilities for a sample size n=12
pik=inclusionprobabilities(a,12)
pik
}
}
