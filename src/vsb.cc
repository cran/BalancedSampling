#include <stddef.h>

#include <Rcpp.h>

#include "KDStoreClass.h"
#include "KDTreeClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.vsb0_cpp)]]
double vsb0_cpp(
  Rcpp::NumericVector &probs,
  Rcpp::NumericVector &ys,
  Rcpp::NumericMatrix &xs,
  size_t treeBucketSize,
  int treeMethod
) {
  size_t N = xs.ncol();
  size_t p = xs.nrow();

  double* xx = REAL(xs);
  double* yp = new double[N];

  KDTree tree(xx, N, p, treeBucketSize, IntToKDTreeSplitMethod(treeMethod));
  KDStore store(N, 1);

  for (size_t i = 0; i < N; i++)
    yp[i] = ys[i] / probs[i];

  double result = 0.0;

  for (size_t i = 0; i < N; i++) {
    tree.FindNeighbours(&store, i);
    size_t len = store.GetSize();

    double localMean = yp[i];

    for (size_t j = 0; j < len; j++)
      localMean += yp[store.neighbours[j]];

    localMean = yp[i] - localMean / (double)(len + 1);
    result += (double)(len + 1) / (double)(len) * (localMean * localMean);
  }

  delete[] yp;

  return result;
}

// [[Rcpp::export(.vsbn_cpp)]]
double vsbn_cpp(
  Rcpp::NumericVector &probs,
  Rcpp::NumericVector &ys,
  Rcpp::NumericMatrix &xs,
  size_t n,
  size_t treeBucketSize,
  int treeMethod
) {
  if (n == (size_t)0)
    throw std::range_error("n must be >= 1");

  size_t N = xs.ncol();
  size_t p = xs.nrow();

  double* xx = REAL(xs);
  double* yp = new double[N];

  KDTree tree(xx, N, p, treeBucketSize, IntToKDTreeSplitMethod(treeMethod));
  KDStore store(N, n);

  for (size_t i = 0; i < N; i++)
    yp[i] = ys[i] / probs[i];

  double result = 0.0;

  for (size_t i = 0; i < N; i++) {
    tree.FindNeighbours(&store, i);
    size_t len = store.GetSize();

    double localMean = yp[i];

    for (size_t j = 0; j < len; j++)
      localMean += yp[store.neighbours[j]];

    localMean = yp[i] - localMean / (double)(len + 1);
    result += (double)(len + 1) / (double)(len) * (localMean * localMean);
  }

  delete[] yp;

  return result;
}


