#include "lin.h"
#include <math.h>

// TODO: replace this comment with a short explanation of this algorithm
_Bool lin_solve(double *m, double *v, int n) {
    for (int i = 0; i < n; i++) {
        // TODO: replace this comment with a one-line explanation of this step
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

        // TODO: replace this comment with a one-line explanation of this step
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

        // TODO: replace this comment with a one-line explanation of this step
        double const diag = m[i*n + i];
        for (int c = i; c < n; c++) {
            m[i*n + c] /= diag;
        }
        {
            v[i] /= diag;
        }

        // TODO: replace this comment with a one-line explanation of this step
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
