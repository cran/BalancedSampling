#include <stddef.h>
#include <stdexcept>

#include <Rcpp.h>

#include "IndexListClass.h"
#include "KDTreeClass.h"
#include "LpmClass.h"

// [[Rcpp::export(.hlpm2_cpp)]]
Rcpp::IntegerMatrix hlpm2_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &sizes,
  size_t treeBucketSize,
  int treeMethod,
  double eps
) {
  size_t sn = sizes.length();
  size_t N = x.ncol();
  size_t p = x.nrow();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob an x does not match");

  Lpm lpm(
    LpmMethod::lpm2,
    REAL(prob),
    REAL(x),
    N,
    p,
    eps,
    treeBucketSize,
    treeMethod
  );

  // Make a copy of the tree for reuse later
  KDTree* orgTree = lpm.tree->Copy();
  IndexList* orgIdx = lpm.idx;

  // Run the algorithm to get a base sample
  lpm.Run();

  orgIdx->Reset();
  delete lpm.tree;

  // Set the return matrix
  size_t orgSampleSize = lpm.sample.size();
  Rcpp::IntegerMatrix sample(orgSampleSize, 2);

  for (size_t i = 0, j = 0; i < N; i++) {
    if (j < orgSampleSize && i == lpm.sample[j] - 1) {
      sample(j, 0) = lpm.sample[j];
      sample(j, 1) = sn;
      j += 1;
    } else {
      orgTree->RemoveUnit(i);
      orgIdx->Erase(i);
    }
  }

  size_t remainingSize = orgSampleSize;

  for (size_t i = 0; i < sn - 1; i++) {
    double subprob = (double)sizes[i] / (double)remainingSize;
    lpm.tree = orgTree->Copy();
    lpm.idx = orgIdx->CopyLen();
    lpm.sample.resize(0);

    for (size_t j = 0; j < orgIdx->Length(); j++) {
      lpm.probabilities[orgIdx->Get(j)] = subprob;
    }

    lpm.Run();

    // Remove all selected unit from orgTree and orgIdx, and set their subsample
    for (size_t j = 0, k = 0; j < orgSampleSize && k < lpm.sample.size(); j++) {
      if ((size_t)sample(j, 0) != lpm.sample[k])
        continue;

      size_t id = lpm.sample[k] - 1;
      orgTree->RemoveUnit(id);
      orgIdx->Erase(id);
      sample(j, 1) = i + 1;
      k += 1;
    }

    remainingSize -= lpm.sample.size();
    delete lpm.tree;
    delete lpm.idx;
  }

  lpm.tree = nullptr;
  lpm.idx = nullptr;
  delete orgTree;
  delete orgIdx;

  return sample;
}

