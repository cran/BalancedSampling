#include <stddef.h>
#include <stdexcept>

#include <Rcpp.h>

#include "CubeStratifiedClass.h"

// [[Rcpp::export(.cube_stratified_cpp)]]
Rcpp::IntegerVector cube_stratified_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &strata,
  const double eps
) {
  size_t N = x.nrow();
  size_t p = x.ncol();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob and x does not match");
  if (N != (size_t)strata.length())
    throw std::range_error("strata and x does not match");

  CubeStratified cube(
    INTEGER(strata),
    REAL(prob),
    REAL(x),
    N,
    p,
    eps
  );

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample_.begin(), cube.sample_.end());

  return sample;
}

// [[Rcpp::export(.lcube_stratified_cpp)]]
Rcpp::IntegerVector lcube_stratified_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &xbalance,
  Rcpp::NumericMatrix &xspread,
  Rcpp::IntegerVector &strata,
  const size_t bucketSize,
  const int method,
  const double eps
) {
  size_t N = xbalance.nrow();
  size_t p = xbalance.ncol();
  size_t pxs = xspread.nrow();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob and x does not match");
  if (N != (size_t)strata.length())
    throw std::range_error("strata and x does not match");
  if (N != (size_t)xspread.ncol())
    throw std::range_error("xspread and xbal does not match");

  CubeStratified cube(
    INTEGER(strata),
    REAL(prob),
    REAL(xbalance),
    N,
    p,
    eps,
    REAL(xspread),
    pxs,
    bucketSize,
    method
  );

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample_.begin(), cube.sample_.end());

  return sample;
}
