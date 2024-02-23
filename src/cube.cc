#include <stddef.h>
#include <stdexcept>

#include <Rcpp.h>

#include "CubeClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.cube_cpp)]]
Rcpp::IntegerVector cube_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  double eps
) {
  size_t N = x.nrow();
  size_t p = x.ncol();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob and x does not match");

  Cube cube(REAL(prob), REAL(x), N, p, eps);

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample.begin(), cube.sample.end());

  return sample;
}

// [[Rcpp::export(.lcube_cpp)]]
Rcpp::IntegerVector lcube_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &xbal,
  Rcpp::NumericMatrix &xspread,
  size_t treeBucketSize,
  int treeMethod,
  double eps
) {
  size_t N = xbal.nrow();
  size_t pbal = xbal.ncol();
  size_t pspread = xspread.nrow();

  if (N != (size_t)xspread.ncol())
    throw std::invalid_argument("xbal and xspread does not match");
  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob and x does not match");

  Cube cube(
    REAL(prob),
    REAL(xbal),
    N,
    pbal,
    eps,
    REAL(xspread),
    pspread,
    treeBucketSize,
    treeMethod
  );

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample.begin(), cube.sample.end());

  return sample;
}
