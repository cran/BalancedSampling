.eps_check = function(eps) {
  if (eps < 0.0 || 1e-4 < eps || !is.numeric(eps) || length(eps) != 1)
    stop("'eps' must be in [0.0, 1e-4]");
}

.kdtree_method_check = function(type, bucketSize) {
  if (type == "kdtree0") {
    method = 0;
  } else if (type == "kdtree1") {
    method = 1;
  } else if (type == "kdtree2") {
    method = 2;
  } else if (type == "notree") {
    method = 0;
  } else {
    stop("'type' must be 'kdtree0', 'kdtree1', 'kdtree2', or 'notree'");
  }

  return(method);
}

.kdtree_bucket_check = function(N, type, bucketSize) {
  if (type == "notree") {
    bucketSize = N;
  } else {
    bucketSize = as.numeric(bucketSize)[1L];
  }

  if (length(bucketSize) != 1 || bucketSize < 1 || bucketSize %% 1 != 0)
    stop("'bucketSize' must be integer > 0");

  return(bucketSize);
}

.prob_check = function(prob, N) {
  if (length(prob) != N)
    stop("the size of 'prob' and 'x' does not match");

  return(as.numeric(prob));
}

.prob_integer_test = function(prob, N) {
  if (length(prob) == 1 && is.numeric(prob)) {
    if (prob < 1 || prob > N || prob %% 1 != 0)
      stop("'prob' must be a vector of probabilities or a single integer in [0, N]");

    return(TRUE);
  }

  return(FALSE);
}

.prob_expand = function(prob, N) {
  if (.prob_integer_test(prob, N)) {
    prob = rep(prob / N, N);
    return(prob);
  }

  return(.prob_check(prob, N));
}

.strata_check = function(strata, N) {
  if (length(strata) != N)
    stop("the size of 'strata' and 'x' does not match");

  return(as.integer(strata))
}
