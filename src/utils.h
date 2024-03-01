#ifndef BSUTILS_HEADER
#define BSUTILS_HEADER

#include <stddef.h>

inline bool Probability1(double p, double eps) {
  return (p) >= 1.0 - (eps);
};
inline bool Probability1(size_t p, size_t N) {
  return (p) == (N);
};

inline bool Probability0(double p, double eps) {
  return (p) <= (eps);
};
inline bool Probability0(size_t p) {
  return (p) == (0);
};

inline bool ProbabilityInt(double p, double eps) {
  return (p) <= (eps) || (p) >= 1.0 - (eps);
};
inline bool ProbabilityInt(size_t p, size_t N) {
  return (p) == (0) || (p) == (N);
};

#endif
