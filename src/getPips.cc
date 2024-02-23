#include <stddef.h>

#include <Rcpp.h>

// [[Rcpp::export(.getpips_cpp)]]
Rcpp::NumericVector getpips_cpp(
  Rcpp::NumericVector &x,
  size_t n
) {
  // Assume < 0 and equality has already been checked
  size_t N = x.length();
  Rcpp::NumericVector prob(N);
  size_t* index = new size_t[N];
  size_t one = 0;
  size_t oneplus = 0;
  size_t oneminus = 0;



  double xsum = 0.0;
  for (size_t i = 0; i < N; i++) {
    if (x[i] < 0.0)
      throw std::range_error("elements in x must be >= 0.0");

    xsum += x[i];
  }

  double temp = ((double) n) / xsum;
  for (size_t i = 0; i < N; i++) {
    prob[i] = x[i] * temp;

    if (prob[i] == 1.0) {
      one += 1;
    } else if (prob[i] > 1.0) {
      prob[i] = 1.0;
      oneplus += 1;
    } else {
      index[oneminus] = i;
      oneminus += 1;
    }
  }

  if (oneplus == 0) {
    delete[] index;
    return prob;
  }

  size_t onesum = one + oneplus;
  while (oneplus > 0) {
    xsum = 0.0;
    for (size_t i = 0; i < oneminus; i++)
      xsum += x[index[i]];

    temp = ((double) (n - onesum)) / xsum;

    size_t toneminus = oneminus;
    one = 0;
    oneplus = 0;
    oneminus = 0;


    for (size_t i = 0; i < toneminus; i++) {
      prob[index[i]] = x[index[i]] * temp;

      if (prob[index[i]] == 1.0) {
        one += 1;
      } else if (prob[index[i]] > 1.0) {
        prob[index[i]] = 1.0;
        oneplus += 1;
      } else {
        index[oneminus] = index[i];
        oneminus += 1;
      }
    }

    onesum += one + oneplus;
  }

  delete[] index;
  return prob;
}
