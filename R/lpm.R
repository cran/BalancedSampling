# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' The (Local) Pivotal Methods
#'
#' @description
#' Selects spatially balanced samples with prescribed inclusion probabilities
#' from a finite population using the Local Pivotal Method 1 (LPM1).
#'
#' @details
#' If `prob` sum to an integer n, a fixed sized sample (n) will be produced.
#' For `spm` and `rpm`, `prob` must be a vector of inclusion probabilities.
#' If equal inclusion probabilities is wanted, this can be produced by
#' `rep(n / N, N)`.
#'
#' The available pivotal methods are:
#' - `lpm1`: The Local Pivotal Mehtod 1 (Grafström et al., 2012).
#' Updates only units which are mutual nearest neighbours.
#' Selects such a pair at random.
#' - `lpm2`, `lpm`: The Local Pivotal Method 2 (Grafström et al., 2012).
#' Selects a unit at random, which competes with this units nearest neighbour.
#' - `lpm1s`: The Local Pivotal Method 1 search: (Prentius, 2023).
#' Updates only units which are mutual nearest neighbours.
#' Selects such a pair by branching the remaining units, giving higher
#' probabilities to update a pair with a long branch.
#' This changes the algorithm of lpm1, but makes it faster.
#' - `spm`: The Sequential Pivotal Method.
#' Selects the two units with smallest indices to compete against each other.
#' If the list is ordered, the algorithm is similar to systematic sampling.
#' - `rpm`: The Random Pivotal Method.
#' Selects two units at random to compete against each other.
#' Produces a design with high entropy.
#'
#'
#' @templateVar xspread x
#' @templateVar integerprob TRUE
#' @template sampling_template
#' @template kdtrees_template
#' @template x_template
#' @template probs_template
#'
#' @references
#' Friedman, J. H., Bentley, J. L., & Finkel, R. A. (1977).
#' An algorithm for finding best matches in logarithmic expected time.
#' ACM Transactions on Mathematical Software (TOMS), 3(3), 209-226.
#'
#' Deville, J.-C., &  Tillé, Y. (1998).
#' Unequal probability sampling without replacement through a splitting method.
#' Biometrika 85, 89-101.
#'
#' Maneewongvatana, S., & Mount, D. M. (1999, December).
#' It’s okay to be skinny, if your friends are fat.
#' In Center for geometric computing 4th annual workshop on computational geometry (Vol. 2, pp. 1-8).
#'
#' Chauvet, G. (2012).
#' On a characterization of ordered pivotal sampling.
#' Bernoulli, 18(4), 1320-1340.
#'
#' Grafström, A., Lundström, N.L.P. & Schelin, L. (2012).
#' Spatially balanced sampling through the Pivotal method.
#' Biometrics 68(2), 514-520.
#'
#' Lisic, J. J., & Cruze, N. B. (2016, June).
#' Local pivotal methods for large surveys.
#' In Proceedings of the Fifth International Conference on Establishment Surveys.
#'
#' Prentius, W. (2023)
#' Manuscript.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = lpm2(prob, x);
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
#'   s = lpm2(prob, x);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#'
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' lpm1(prob, x);
#' lpm2(prob, x);
#' lpm1s(prob, x);
#' spm(prob);
#' rpm(prob);
#' }
lpm = function(
  prob,
  x,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  lpm2(prob, x, type, bucketSize, eps);
}

#' Local Pivotal Method 1 (LPM1)
#'
#' @describeIn lpm
#'
lpm1 = function(
  prob,
  x,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = t(as.matrix(x));
  } else {
    x = t(x);
  }

  N = dim(x)[2L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  .eps_check(eps);

  if (.prob_integer_test(prob, N)) {
    result = .lpm_int_cpp(1, prob, x, bucketSize, method);
  } else {
    prob = .prob_check(prob, N);
    result = .lpm_cpp(1, prob, x, bucketSize, method, eps);
  }

  return(result);
}

#' Local Pivotal Method 2 (LPM2)
#'
#' @describeIn lpm
#'
lpm2 = function(
  prob,
  x,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = t(as.matrix(x));
  } else {
    x = t(x);
  }

  N = dim(x)[2L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  .eps_check(eps);

  if (.prob_integer_test(prob, N)) {
    result = .lpm_int_cpp(2, prob, x, bucketSize, method);
  } else {
    prob = .prob_check(prob, N);
    result = .lpm_cpp(2, prob, x, bucketSize, method, eps);
  }

  return(result);
}

#' Local Pivotal Method 1 search (LPM1s)
#'
#' @describeIn lpm
#'
lpm1s = function(
  prob,
  x,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = t(as.matrix(x));
  } else {
    x = t(x);
  }

  N = dim(x)[2L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  .eps_check(eps);

  if (.prob_integer_test(prob, N)) {
    result = .lpm_int_cpp(3, prob, x, bucketSize, method);
  } else {
    prob = .prob_check(prob, N);
    result = .lpm_cpp(3, prob, x, bucketSize, method, eps);
  }

  return(result);
}

#' Sequential Pivotal Method (SPM)
#'
#' @describeIn lpm
#'
spm = function(
  prob,
  eps = 1e-12
) {
  if (length(prob) == 1)
    stop("'prob' must be a vector of probabilities");

  prob = as.numeric(prob);
  .eps_check(eps);

  result = .spm_cpp(prob, eps);

  return(result);
}

#' Random Pivotal Method (RPM)
#'
#' @describeIn lpm
#'
rpm = function(
  prob,
  eps = 1e-12
) {
  if (length(prob) == 1)
    stop("'prob' must be a vector of probabilities");

  prob = as.numeric(prob);
  .eps_check(eps);

  result = .rpm_cpp(prob, eps);

  return(result);
}
