#include <cmath>
#include <stddef.h>

#include <RcppArmadillo.h>

#include "KDStoreClass.h"
#include "KDTreeClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.sb_localbalance_cpp)]]
double sb_localbalance_cpp(
  arma::vec &prob,
  arma::mat &x,
  arma::uvec &sample,
  size_t treeBucketSize,
  size_t treeMethod
) {
  size_t N = x.n_cols;
  size_t p = x.n_rows;
  size_t n = sample.n_elem;
  size_t p1 = p + 1;

  if (prob.n_elem != N)
    throw std::out_of_range("size of 'prob' does not match x");

  // xd contains the diffs between sample units (prob weighted) and voronoi set
  arma::mat xd(p1, n, arma::fill::none);
  arma::mat xs(p, n, arma::fill::none);
  arma::mat qq(p1, p1, arma::fill::zeros);

  for (size_t i = 0; i < n; i++) {
    if (sample[i] < (size_t)1 || sample[i] > (size_t)N)
      throw std::range_error("'sample' must contain unit indices");

    size_t idx = sample[i] - 1;
    double probFactor = (1.0 - prob(idx)) / prob(idx);

    xs.col(i) = x.col(idx);

    xd.col(i).tail(p) = x.col(idx) * probFactor;
    xd(0, i) = probFactor;
  }

  KDTree tree(xs.memptr(), n, p, treeBucketSize, IntToKDTreeSplitMethod(treeMethod));
  KDStore store(n, 1);

  // Prepare xd
  for (size_t i = 0; i < N; i++) {
    arma::vec unit(p1, arma::fill::none);
    unit(0) = 1.0;
    unit.tail(p) = x.col(i);

    // Add to the qq-matrix
    qq += unit * unit.t();

    // Check if unit is in sample
    size_t break_outer = 0;

    for (size_t j = 0; j < n; j++) {
      if (i == sample[j] - 1) {
        break_outer = 1;
        break;
      }
    }

    if (break_outer == 1) continue;

    // Find neighbours to population unit in sample
    double* unitptr = x.colptr(i);
    tree.FindNeighbours(&store, unitptr);
    size_t storeSize = store.GetSize();

    // If only one neighbour, all is fine and dandy
    if (storeSize == 1) {
      size_t sampleUnit = store.neighbours[0];
      xd.col(sampleUnit) -= unit;
      continue;
    }

    // Handle tied neighbours, by evenly splitting the unit amongst the closest
    unit /= double(storeSize);
    for (size_t j = 0; j < storeSize; j++) {
      size_t sampleUnit = store.neighbours[j];
      xd.col(sampleUnit) -= unit;
    }
  }

  // Change qq to inverse
  qq = qq.i();

  // Sum the results
  double result = 0.0;

  for (size_t i = 0; i < n; i++) {
    arma::vec temp = xd.col(i);
    result += arma::as_scalar(temp.t() * qq * temp);
  }

  return std::sqrt(result / double(N));
}

