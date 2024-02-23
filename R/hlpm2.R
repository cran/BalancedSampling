# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Hierarchical Local Pivotal Method 2
#'
#' @description
#' Selects an initial sample using the [lpm2()], and then splits this sample into
#' subsamples of given `sizes` using successive, hierarchical selection with
#' the [lpm2()].
#' The method is used to select several subsamples, such that each subsample, and
#' the combination (i.e. the union of all subsamples), is spatially balanced.
#'
#' @details
#' The inclusion probabilities `prob` _must_ sum to an integer n.
#' The sizes of the subsamples `sum(sizes)` _must_ sum to the same integer n.
#'
#' @templateVar xspread x
#' @template sampling_template
#' @template kdtrees_template
#' @template x_template
#' @template probs_template
#'
#' @param sizes A vector of integers containing the sizes of the subsamples.
#' `sum(sizes) = sum(prob)` must hold.
#'
#' @return A matrix with the population indices of the combined sample in the
#' first column, and the associated subsample in the second column.
#'
#' @references
#' Friedman, J. H., Bentley, J. L., & Finkel, R. A. (1977).
#' An algorithm for finding best matches in logarithmic expected time.
#' ACM Transactions on Mathematical Software (TOMS), 3(3), 209-226.
#'
#' Maneewongvatana, S., & Mount, D. M. (1999, December).
#' It’s okay to be skinny, if your friends are fat.
#' In Center for geometric computing 4th annual workshop on computational geometry (Vol. 2, pp. 1-8).
#'
#' Grafström, A., Lundström, N.L.P. & Schelin, L. (2012).
#' Spatially balanced sampling through the Pivotal method.
#' Biometrics 68(2), 514-520.
#'
#' Lisic, J. J., & Cruze, N. B. (2016, June).
#' Local pivotal methods for large surveys.
#' In Proceedings of the Fifth International Conference on Establishment Surveys.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' sizes = c(10, 20, 30, 40);
#' s = hlpm2(prob, x, sizes);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#' }
#'
hlpm2 = function(
  prob,
  x,
  sizes,
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
  prob = .prob_check(prob, N);

  probsum = sum(prob);

  if (probsum %% 1 != 0)
    stop("'prob' must sum to an integer");

  if (probsum != sum(sizes))
    stop("'sizes' must sum to an integer same as the sum of 'prob'");

  result = .hlpm2_cpp(prob, x, sizes, bucketSize, method, eps);

  return(result);
}
