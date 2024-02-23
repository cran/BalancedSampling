#include <stddef.h>

#include "ReducedRowEchelonForm.h"
#include "utils.h"

void ReducedRowEchelonForm(
  double *mat,
  const size_t rowCount,
  const size_t colCount
) {
  size_t lead = 0;

  for (size_t r = 0; r < rowCount; r++) {
    if (colCount <= lead)
      return;

    size_t i = r;

    while (mat[MatrixIdx(i, lead, colCount)] == 0) {
      i += 1;

      if (i == rowCount) {
        i = r;
        lead += 1;

        if (colCount == lead)
          return;
      }
    }

    double *br = mat + MatrixIdx(r, (size_t)0, colCount);

    if (i != r) {
      double *bi = mat + MatrixIdx(i, (size_t)0, colCount);

      // Swap rows
      for (size_t k = 0; k < colCount; k++) {
        double temp = bi[k];
        bi[k] = br[k];
        br[k] = temp;
      }
    }

    // Divide by lead, assuming all is 0 before lead
    /*if (br[lead] != 0.0)*/ {
      double temp = br[lead];
      br[lead] = 1.0;
      for (size_t k = lead + 1; k < colCount; k++) {
        br[k] /= temp;
      }
    }

    // Remove row r from all outher rows
    for (size_t j = 0; j < rowCount; j++) {
      if (j == r)
        continue;

      double temp = mat[MatrixIdx(j, lead, colCount)];
      mat[MatrixIdx(j, lead, colCount)] = 0.0;
      for (size_t k = lead + 1; k < colCount; k++) {
        mat[MatrixIdx(j, k, colCount)] -= br[k] * temp;
      }
    }

    lead += 1;
  }

  return;
}
