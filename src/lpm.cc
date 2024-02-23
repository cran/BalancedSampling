#include <stddef.h>
#include <stdexcept>

#include <Rcpp.h>

#include "LpmClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.lpm_cpp)]]
Rcpp::IntegerVector lpm_cpp(
  int lpMethod,
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

  Lpm lpm(
    IntToLpmMethod(lpMethod),
    REAL(prob),
    REAL(x),
    N,
    p,
    eps,
    treeBucketSize,
    treeMethod
  );

  lpm.Run();

  Rcpp::IntegerVector sample(lpm.sample.begin(), lpm.sample.end());

  return sample;
}

// [[Rcpp::export(.lpm_int_cpp)]]
Rcpp::IntegerVector lpm_int_cpp(
  int lpMethod,
  size_t n,
  Rcpp::NumericMatrix &x,
  size_t treeBucketSize,
  int treeMethod
) {
  size_t N = x.ncol();
  size_t p = x.nrow();

  Lpm lpm(
    IntToLpmMethod(lpMethod),
    n,
    REAL(x),
    N,
    p,
    treeBucketSize,
    treeMethod
  );

  lpm.Run();

  Rcpp::IntegerVector sample(lpm.sample.begin(), lpm.sample.end());

  return sample;
}

// [[Rcpp::export(.rpm_cpp)]]
Rcpp::IntegerVector rpm_cpp(
  Rcpp::NumericVector &prob,
  double eps
) {
  size_t N = prob.length();

  Lpm lpm(
    LpmMethod::rpm,
    REAL(prob),
    nullptr,
    N,
    0,
    eps,
    40,
    2
  );

  lpm.Run();

  Rcpp::IntegerVector sample(lpm.sample.begin(), lpm.sample.end());

  return sample;
}

// [[Rcpp::export(.spm_cpp)]]
Rcpp::IntegerVector spm_cpp(
  Rcpp::NumericVector &prob,
  double eps
) {
  size_t N = prob.length();

  Lpm lpm(
    LpmMethod::spm,
    REAL(prob),
    nullptr,
    N,
    0,
    eps,
    40,
    2
  );

  lpm.Run();

  Rcpp::IntegerVector sample(lpm.sample.begin(), lpm.sample.end());

  return sample;
}
