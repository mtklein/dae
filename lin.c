#include "lin.h"
#include <math.h>

// Solve a linear system using Gauss-Jordan elimination with partial pivoting.
_Bool lin_solve(double *m, double *v, int n) {
    for (int i = 0; i < n; i++) {
        // Locate the row with the largest coefficient for this column.
        double max = fabs(m[i*n + i]);
        int    piv = i;
        for (int r = i+1; r < n; r++) {
            double const val = fabs(m[r*n + i]);
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
                double const swap = m[  i*n + c];
                m[  i*n + c]      = m[piv*n + c];
                m[piv*n + c]      = swap;
            }
            {
                double const swap = v[  i];
                v[  i]            = v[piv];
                v[piv]            = swap;
            }
        }

        // Scale the pivot row so the pivot column becomes 1.
        double const diag = m[i*n + i];
        for (int c = i; c < n; c++) {
            m[i*n + c] /= diag;
        }
        {
            v[i] /= diag;
        }

        // Eliminate the pivot column from all other rows.
        for (int r = 0; r < n; r++) {
            if (r != i) {
                double const factor = m[r*n + i];
                for (int c = i; c < n; c++) {
                    m[r*n + c] -= factor * m[i*n + c];
                }
                v[r] -= factor * v[i];
            }
        }
    }
    return (_Bool)1;
}
