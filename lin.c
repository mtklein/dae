#include "lin.h"
#include <math.h>

// Solve a linear system using Gauss-Jordan elimination with partial pivoting.
_Bool lin_solve(double *M, double *V, int n) {
    for (int i = 0; i < n; i++) {
        // Locate the row with the largest coefficient for this column.
        double max = fabs(M[i*n + i]);
        int    piv = i;
        for (int r = i+1; r < n; r++) {
            double const val = fabs(M[r*n + i]);
            if (max < val) {
                max = val;
                piv = r;
            }
        }
        if (max == 0.0) {
            return (_Bool)0;
        }

        // Swap the current row with the pivot row in the matrix and vector.
        if (piv != i) {
            for (int c = i; c < n; c++) {
                double const swap = M[  i*n + c];
                M[  i*n + c]      = M[piv*n + c];
                M[piv*n + c]      = swap;
            }
            {
                double const swap = V[  i];
                V[  i]            = V[piv];
                V[piv]            = swap;
            }
        }

        // Scale the pivot row so the pivot column becomes 1.
        double const diag = M[i*n + i];
        for (int c = i; c < n; c++) {
            M[i*n + c] /= diag;
        }
        {
            V[i] /= diag;
        }

        // Eliminate the pivot column from all other rows.
        for (int r = 0; r < n; r++) {
            if (r != i) {
                double const factor = M[r*n + i];
                for (int c = i; c < n; c++) {
                    M[r*n + c] -= factor * M[i*n + c];
                }
                V[r] -= factor * V[i];
            }
        }
    }
    return (_Bool)1;
}
