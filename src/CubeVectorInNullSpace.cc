#include <stddef.h>
#include <stdexcept>

#include "CubeVectorInNullSpace.h"
#include "utils-matrix.h"

// res_vec assumed to have length (ncolumns)
// mat Matrix assumed to be in Row Major form, dimensions (ncolumns-1) * (ncolumns)
void CubeVectorInNullSpace(
  double *res_vec,
  double *mat,
  const size_t ncolumns
) {
  size_t nrows = ncolumns - 1;
  if (ncolumns < 2) {
    throw std::range_error("nrows and ncolumns must be >= 2");
  }
  if (mat[MatrixIdxRM((size_t)0, (size_t)0, ncolumns)] == 0.0) {
    throw std::range_error("no null basis exists");
  }

  // If the obs. are linearly independent, we can take the fast route
  if (mat[MatrixIdxRM(nrows-1, nrows-1, ncolumns)] == 1.0) {
    res_vec[ncolumns-1] = 1.0;

    for (size_t i = 0; i < nrows; i++) {
      res_vec[i] = -mat[MatrixIdxRM(i, ncolumns-1, ncolumns)];
    }

    return;
  }

  // If any obs is linearly dependent, we must take the slow route
  // We can assume that the first entry in mat is 1.0
  for (size_t k = 1; k < ncolumns; k++) res_vec[k] = (k % 2 == 0) ? -1.0 : 1.0;

  for (size_t i = 0; i < nrows; i++) {
    size_t lead = 0;

    for (; lead < ncolumns; lead++) {
      if (mat[MatrixIdxRM(i, lead, ncolumns)] == 1.0) break;
    }

    if (lead >= ncolumns) continue;

    res_vec[lead] = 0.0;
    for (size_t k = lead + 1; k < ncolumns; k++) {
      res_vec[lead] -= res_vec[k] * mat[MatrixIdxRM(i, k, ncolumns)];
    }
  }
}

