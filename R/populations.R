.pop_poisson = function(parents, meanChildren, distribution) {
  NP = nrow(parents);
  children = rpois(NP, meanChildren);
  ccs = cumsum(children);

  # Draw the population
  population = matrix(distribution(sum(children)), sum(children), ncol(parents), TRUE);

  # Center the children around their parents
  start = 1L;
  stop = 0L;

  for (i in seq_len(NP)) {
    stop = ccs[i];

    population[start:stop, ] = sweep(
      population[start:stop, ],
      2,
      parents[i, ],
      '+'
    );

    start = stop + 1L;
  }

  return(population);
}

.pop_mirror = function(population, minv, maxv) {
  N = nrow(population);
  P = ncol(population);

  if (length(minv) == 1)
    minv = rep(minv, P)
  if (length(maxv) == 1)
    maxv = rep(maxv, P)

  diff = maxv - minv;
  temp = 0.0;
  flor = 0;

  for (i in seq_len(N)) {
    for (j in seq_len(P)) {
      if (population[i, j] >= minv[j] && population[i, j] <= maxv[j]) {
        next;
      }

      temp = population[i, j] - minv[j];
      flor = floor(temp / diff[j]);

      if (flor %% 2 == 0) {
        population[i, j] = temp %% diff[j] + minv[j];
      } else {
        population[i, j] = diff[j] - (temp %% diff[j]) + minv[j];
      }
    }
  }

  return(population);
}

#' Generate populations
#'
#' @importFrom stats rnorm rpois runif
#' @description
#' Generate uniform and poisson cluster process populations
#'
#' If `from` and `to` is used with `genpopPoisson` together with `mirror`, the
#' population will be bounded within these values.
#' For the `genpopUniform`, these numbers represent the minimum and maximum
#' values of the uniform distribution.
#'
#' @param size The size of the population
#' @param dims The number of auxiliary variables
#' @param from A number or a vector of size `dims` with the minimum values
#' @param to A number or a vector of size `dims` with the maximum values
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' x = genpopUniform(120, 2L);
#' N = nrow(x);
#' n = 60;
#' prob = rep(n / N, N);
#' s = lpm2(prob, x);
#' b = sb(prob, x, s);
#' }
#'
genpopUniform = function(size, dims = 2L, from = 0.0, to = 1.0) {
  population = matrix(
    runif(size * dims, from, to),
    size,
    dims,
    TRUE
  );

  return(population);
}

#' @describeIn genpopUniform Poisson cluster process
#' @param parents The number of parent locations
#' @param children A number or a vector of size `parents` with the mean number of
#' children to be spawned.
#' @param distribution A function taking a number as a variable, returning the
#' offset from the parent location.
#' @param mirror If `TRUE`, the population is mirrored to be inside `from` and `to`.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' x = genpopPoisson(70, 50, 2L);
#' N = nrow(x);
#' n = 60;
#' prob = rep(n / N, N);
#' s = lpm2(prob, x);
#' b = sb(prob, x, s);
#' }
#'
genpopPoisson = function(
  parents,
  children,
  dims = 2L,
  from = 0.0,
  to = 1.0,
  distribution = function(n) rnorm(n, 0.0, 0.02),
  mirror = TRUE
) {
  parpop = genpopUniform(parents, dims, from, to);
  population = .pop_poisson(parpop, children, distribution);

  if (mirror == TRUE) {
    population = .pop_mirror(population, from, to);
  }

  return(population);
}
