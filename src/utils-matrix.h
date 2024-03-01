#ifndef BSUTILS_MATRIX_HEADER
#define BSUTILS_MATRIX_HEADER

#include <stddef.h>

// ROW MAJOR
inline int MatrixIdxRM(int row, int col, int pcols) {
  return row * pcols + col;
}
inline size_t MatrixIdxRM(size_t row, size_t col, size_t pcols) {
  return row * pcols + col;
}

// COL MAJOR
inline int MatrixIdxCM(int row, int col, int nrows) {
  return col * nrows + row;
}
inline size_t MatrixIdxCM(size_t row, size_t col, size_t nrows) {
  return col * nrows + row;
}

#endif

