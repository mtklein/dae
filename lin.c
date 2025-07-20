#include "lin.h"
#include <math.h>

int lin_solve(double *matrix, double *vec, int n) {
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        double max = fabs(matrix[i * n + i]);
        for (int r = i + 1; r < n; ++r) {
            double v = fabs(matrix[r * n + i]);
            if (v > max) {
                max = v;
                pivot = r;
            }
        }
        if (max == 0.0) {
            return -1;
        }
        if (pivot != i) {
            for (int c = i; c < n; ++c) {
                double tmp = matrix[i * n + c];
                matrix[i * n + c] = matrix[pivot * n + c];
                matrix[pivot * n + c] = tmp;
            }
            double tmp = vec[i];
            vec[i] = vec[pivot];
            vec[pivot] = tmp;
        }
        double diag = matrix[i * n + i];
        for (int c = i; c < n; ++c) {
            matrix[i * n + c] /= diag;
        }
        vec[i] /= diag;
        for (int r = 0; r < n; ++r) {
            if (r == i) {
                continue;
            }
            double factor = matrix[r * n + i];
            for (int c = i; c < n; ++c) {
                matrix[r * n + c] -= factor * matrix[i * n + c];
            }
            vec[r] -= factor * vec[i];
        }
    }
    return 0;
}
