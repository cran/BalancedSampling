# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' The Cube method
#'
#' @description
#' Selects balanced samples with prescribed inclusion probabilities
#' from a finite population using the fast flight Cube Method.
#'
#' @details
#' If `prob` sum to an integer n, and `prob` is included as the first
#' balancing variable, a fixed sized sample (n) will be produced.
#'
#' ## Stratified cube
#' For `cubestratified`, `prob` is automatically inserted as a balancing variable.
#'
#' The stratified version uses the fast flight Cube method and pooling of
#' landing phases.
#'
#' @templateVar xbal x
#' @template sampling_template
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
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = cube(prob, x);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' strata = c(rep(1L, 100), rep(2L, 200), rep(3L, 300), rep(4L, 400));
#' s = cubestratified(prob, x, strata);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' prob = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
#' N = length(prob);
#' x = matrix(runif(N * 2), ncol = 2);
#' ep = rep(0L, N);
#' r = 10000L;
#' for (i in seq_len(r)) {
#'   s = cube(prob, cbind(prob, x));
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
cube = function(
  prob,
  x,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = as.matrix(x);
  }

  N = dim(x)[1L];
  .eps_check(eps);
  prob = .prob_check(prob, N);

  result = .cube_cpp(prob, x, eps);

  return(result);
}

#' Stratified Cube method
#'
#' @describeIn cube
#'
cubestratified = function(
  prob,
  x,
  integerStrata,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = as.matrix(x);
  }

  N = dim(x)[1L];
  .eps_check(eps);
  prob = .prob_check(prob, N);
  strata = .strata_check(integerStrata, N);

  result = .cube_stratified_cpp(prob, x, strata, eps);

  return(result);
}
