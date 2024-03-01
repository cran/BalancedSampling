#include <stddef.h>

#include "ReducedRowEchelonForm.h"
#include "utils-matrix.h"

void ReducedRowEchelonForm(
  double *mat,
  const size_t nrows,
  const size_t ncolumns
) {
  size_t lead = 0;

  for (size_t r = 0; r < nrows; r++) {
    if (ncolumns <= lead)
      return;

    size_t i = r;

    while (mat[MatrixIdxRM(i, lead, ncolumns)] == 0.0) {
      i += 1;

      if (i == nrows) {
        i = r;
        lead += 1;

        if (ncolumns == lead)
          return;
      }
    }

    double *br = mat + MatrixIdxRM(r, (size_t)0, ncolumns);

    if (i != r) {
      double *bi = mat + MatrixIdxRM(i, (size_t)0, ncolumns);

      // Swap rows
      for (size_t k = 0; k < ncolumns; k++) {
        double temp = bi[k];
        bi[k] = br[k];
        br[k] = temp;
      }
    }

    // Divide by lead, assuming all is 0 before lead
    /*if (br[lead] != 0.0)*/ {
      double temp = br[lead];
      br[lead] = 1.0;
      for (size_t k = lead + 1; k < ncolumns; k++) {
        br[k] /= temp;
      }
    }

    // Remove row r from all outher rows
    for (size_t j = 0; j < nrows; j++) {
      if (j == r)
        continue;

      double temp = mat[MatrixIdxRM(j, lead, ncolumns)];
      mat[MatrixIdxRM(j, lead, ncolumns)] = 0.0;
      for (size_t k = lead + 1; k < ncolumns; k++) {
        mat[MatrixIdxRM(j, k, ncolumns)] -= br[k] * temp;
      }
    }

    lead += 1;
  }

  return;
}
