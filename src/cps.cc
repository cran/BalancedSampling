#include <stddef.h>
#include <stdexcept>

#include <Rcpp.h>

#include "CpsClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.cps_cpp)]]
Rcpp::IntegerVector cps_cpp(
  int cpsMethod,
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  size_t treeBucketSize,
  int treeMethod,
  double eps
) {
  size_t N = x.ncol();
  size_t p = x.nrow();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob an x does not match");

  Cps cps(
    IntToCpsMethod(cpsMethod),
    REAL(prob),
    REAL(x),
    N,
    p,
    eps,
    treeBucketSize,
    treeMethod
  );

  cps.Run();

  Rcpp::IntegerVector sample(cps.sample.begin(), cps.sample.end());

  return sample;
}

// [[Rcpp::export(.cps_random_cpp)]]
Rcpp::IntegerVector cps_random_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::NumericVector &random,
  size_t treeBucketSize,
  int treeMethod,
  double eps
) {
  size_t N = x.ncol();
  size_t p = x.nrow();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob an x does not match");
  if (N != (size_t)random.length())
    throw std::invalid_argument("random an x does not match");

  Cps cps(
    CpsMethod::scpscoord,
    REAL(prob),
    REAL(x),
    N,
    p,
    eps,
    treeBucketSize,
    treeMethod
  );

  cps.SetRandomArr(REAL(random));

  cps.Run();

  Rcpp::IntegerVector sample(cps.sample.begin(), cps.sample.end());

  return sample;
}
