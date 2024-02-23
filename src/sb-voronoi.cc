#include <algorithm>
#include <stddef.h>

#include <Rcpp.h>

#include "KDStoreClass.h"
#include "KDTreeClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

//[[Rcpp::export(.sb_voronoi_cpp)]]
double sb_voronoi_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &sample,
  size_t treeBucketSize,
  int treeMethod
) {
  size_t N = x.ncol();
  size_t p = x.nrow();
  size_t n = sample.length();
  double* xx = REAL(x);
  double* xs = new double[n * p];
  double* incl = new double[n];

  for (size_t i = 0; i < n; i++) {
    if (sample[i] < 1 || sample[i] > (int)N)
      throw std::range_error("'sample' must contain unit indices");

    std::copy_n(xx + (sample[i] - 1) * p, p, xs + i * p);
    incl[i] = 0.0;
  }

  KDTree tree(xs, n, p, treeBucketSize, IntToKDTreeSplitMethod(treeMethod));
  KDStore store(n, 1);

  for (size_t i = 0; i < N; i++) {
    double* unit = xx + i * p;
    tree.FindNeighbours(&store, unit);
    size_t len = store.GetSize();

    double ppart = len == (size_t)1 ? prob[i] : prob[i] / (double)len;
    for (size_t j = 0; j < len; j++)
      incl[store.neighbours[j]] += ppart;
  }

  double result = 0.0;
  for (size_t i = 0; i < n; i++) {
    double temp = incl[i] - 1.0;
    result += temp * temp;
  }

  delete[] xs;
  delete[] incl;

  return result / (double)n;
}
