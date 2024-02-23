#ifndef STDUNIFORM_HEADER
#define STDUNIFORM_HEADER

#include <Rcpp.h>

inline double stduniform() {
  double u;
  do {u = R::unif_rand();} while (u < 0.0 || u >= 1.0);
  return u;
};

inline double stduniform(double v) {
  return (stduniform() * v);
};

inline int intuniform(int N) {
  return ((N == 0 || N == 1) ? 0 : (int)stduniform((double)N));
};

inline size_t sizeuniform(size_t N) {
  return ((N == 0 || N == 1) ? 0 : (size_t)stduniform((double)N));
};

/* double stduniform(); */
/* int intuniform(int); */
/* int intuniform(double); */

#endif
