# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Inclusion probabilities proportional-to-size
#'
#' @family utils
#'
#' @description
#' Computes the first-order inclusion probabilties from a vector of positive numbers,
#' for a probabilitiy proportional-to-size design.
#'
#' @param x A vector of positive numbers
#' @param n The wanted sample size
#'
#' @return A vector of inclusion probabilities proportional-to-size
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' x = matrix(runif(N * 2), ncol = 2);
#' prob = getPips(x[, 1], n);
#' s = lpm2(prob, x);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#' }
#'
getPips = function(
  x,
  n
) {
  x = as.numeric(x);
  N = length(x);

  if (length(n) != 1 || n > N || n < 0 || n %% 1 != 0)
    stop("'n' must be integer in [0, N]");

  if (n == N)
    return(rep(1.0, N));
  if (n == 0)
    return(rep(0.0, N));

  result = .getpips_cpp(x, n);

  return(result);
}
