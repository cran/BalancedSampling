# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Variance estimator for spatially balanced samples
#'
#' @family measure
#'
#' @description
#' Variance estimator of HT estimator of population total.
#'
#' @details
#' If `k = 0L`, the variance estimate is constructed by using all units that
#' have the minimum distance.
#'
#' If `k > 0L`, the variance estimate is constructed by using the `k` closest
#' units. If multiple units are located on the border, all are used.
#'
#' @template kdtrees_template
#'
#' @param probs A vector of length n with inclusion probabilities.
#' @param ys A vector of length n containing the target variable.
#' @param xs An n by p matrix of (standardized) auxiliary variables.
#' @param k The number of neighbours to construct the means around.
#'
#' @return The variance estimate.
#'
#' @references
#' Grafström, A., & Schelin, L. (2014).
#' How to select representative samples.
#' Scandinavian Journal of Statistics, 41(2), 277-290.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' y = runif(N);
#' s = lpm2(prob, x);
#' vsb(prob[s], y[s], x[s, ]);
#' vsb(prob[s], y[s], x[s, ], 0L);
#' }
#'

vsb = function(
  probs,
  ys,
  xs,
  k = 3L,
  type = "kdtree2",
  bucketSize = 40
) {
  if (!is.matrix(xs)) {
    xs = t(as.matrix(xs));
  } else {
    xs = t(xs);
  }

  N = dim(xs)[2L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  probs = .prob_check(probs, N);

  if (length(ys) != N)
    stop("the size of 'ys' and 'xs' does not match");

  ys = as.numeric(ys);

  if (k == 0L) {
    result = .vsb0_cpp(probs, ys, xs, bucketSize, method);
  } else {
    result = .vsbn_cpp(probs, ys, xs, k, bucketSize, method);
  }

  return(result);
}
