# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' The Local Cube method
#'
#' @description
#' Selects doubly balanced samples with prescribed inclusion probabilities
#' from a finite population using the Local Cube method.
#'
#' @details
#' If \code{prob} sum to an integer n, and \code{prob} is included as the first
#' balancing variable, a fixed sized sample (n) will be produced.
#'
#' ## Stratified lcube
#' For `lcubestratified`, `prob` is automatically inserted as a balancing variable.
#'
#' The stratified version uses the fast flight Cube method and pooling of
#' landing phases.
#'
#' @templateVar xbal Xbal
#' @templateVar xspread Xspread
#' @template sampling_template
#' @template kdtrees_template
#' @template x_template
#' @template probs_template
#'
#' @param integerStrata An integer vector of length N with stratum numbers.
#'
#' @references
#' Deville, J. C. and Tillé, Y. (2004).
#' Efficient balanced sampling: the cube method.
#' Biometrika, 91(4), 893-912.
#'
#' Chauvet, G. and Tillé, Y. (2006).
#' A fast algorithm for balanced sampling.
#' Computational Statistics, 21(1), 53-62.
#'
#' Chauvet, G. (2009).
#' Stratified balanced sampling.
#' Survey Methodology, 35, 115-119.
#'
#' Grafström, A. and Tillé, Y. (2013).
#' Doubly balanced spatial sampling with spreading and restitution of auxiliary totals.
#' Environmetrics, 24(2), 120-131
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' xspr = matrix(runif(N * 2), ncol = 2);
#' s = lcube(prob, xspr, cbind(prob, x));
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' xspr = matrix(runif(N * 2), ncol = 2);
#' strata = c(rep(1L, 100), rep(2L, 200), rep(3L, 300), rep(4L, 400));
#' s = lcubestratified(prob, xspr, x, strata);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' prob = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
#' N = length(prob);
#' x = matrix(runif(N * 2), ncol = 2);
#' xspr = matrix(runif(N * 2), ncol = 2);
#' ep = rep(0L, N);
#' r = 10000L;
#' for (i in seq_len(r)) {
#'   s = lcube(prob, xspr, cbind(prob, x));
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
lcube = function(
  prob,
  Xspread,
  Xbal,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(Xbal)) {
    Xbal = as.matrix(Xbal);
  }

  if (!is.matrix(Xspread)) {
    Xspread = t(as.matrix(Xspread));
  } else {
    Xspread = t(Xspread);
  }

  N = dim(Xbal)[1L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  .eps_check(eps);
  prob = .prob_check(prob, N);

  if (N != dim(Xspread)[2L])
    stop("the size of 'Xbal' and 'Xspread' does not match");

  result = .lcube_cpp(prob, Xbal, Xspread, bucketSize, method, eps);

  return(result);
}

#' Stratified Local Cube method
#'
#' @describeIn lcube
#'
lcubestratified = function(
  prob,
  Xspread,
  Xbal,
  integerStrata,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(Xbal)) {
    Xbal = as.matrix(Xbal);
  }

  if (!is.matrix(Xspread)) {
    Xspread = t(as.matrix(Xspread));
  } else {
    Xspread = t(Xspread);
  }

  N = dim(Xbal)[1L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  .eps_check(eps);
  prob = .prob_check(prob, N);
  strata = .strata_check(integerStrata, N);

  if (N != dim(Xspread)[2L])
    stop("the size of 'Xbal' and 'Xspread' does not match");

  result = .lcube_stratified_cpp(prob, Xbal, Xspread, strata, bucketSize, method, eps);

  return(result);
}
