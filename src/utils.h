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

inline int MatrixIdx(int row, int col, int pcols) {
  return (row) * (pcols) + (col);
}
inline size_t MatrixIdx(size_t row, size_t col, size_t pcols) {
  return (row) * (pcols) + (col);
}
inline int MatrixIdxRow(int row, int col, int pcols) {
  return (row) * (pcols) + (col);
}
inline size_t MatrixIdxRow(size_t row, size_t col, size_t pcols) {
  return (row) * (pcols) + (col);
}
inline int MatrixIdxCol(int row, int col, int nrows) {
  return (col) * (nrows) + (row);
}
inline size_t MatrixIdxCol(size_t row, size_t col, size_t nrows) {
  return (col) * (nrows) + (row);
}

#endif
