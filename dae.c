#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdlib.h>

int dae_step_euler(struct dae_system const *sys, double step,
                   double *y, double *z, double t) {
    int ny = sys->dim_y;
    int nz = sys->dim_z;

    double *dydt = malloc((size_t)ny * sizeof *dydt);
    if (!dydt) {
        return -1;
    }
    sys->f(y, z, t, dydt, sys->ctx);

    double *y_new = malloc((size_t)ny * sizeof *y_new);
    if (!y_new) {
        free(dydt);
        return -1;
    }
    for (int i = 0; i < ny; ++i) {
        y_new[i] = y[i] + step * dydt[i];
    }

    double *z_guess = malloc((size_t)nz * sizeof *z_guess);
    if (!z_guess) {
        free(y_new);
        free(dydt);
        return -1;
    }
    for (int i = 0; i < nz; ++i) {
        z_guess[i] = z[i];
    }

    for (int iter = 0; iter < 10; ++iter) {
        double *gvec = malloc((size_t)nz * sizeof *gvec);
        if (!gvec) {
            break;
        }
        sys->g(y_new, z_guess, t + step, gvec, sys->ctx);
        double norm = 0.0;
        for (int i = 0; i < nz; ++i) {
            norm += gvec[i] * gvec[i];
        }
        if (norm < 1e-16) {
            free(gvec);
            for (int i = 0; i < ny; ++i) {
                y[i] = y_new[i];
            }
            for (int i = 0; i < nz; ++i) {
                z[i] = z_guess[i];
            }
            free(z_guess);
            free(y_new);
            free(dydt);
            return 0;
        }
        size_t jm = (size_t)nz * (size_t)nz;
        double *J = malloc(jm * sizeof *J);
        double *delta = malloc((size_t)nz * sizeof *delta);
        if (!J || !delta) {
            free(delta);
            free(J);
            free(gvec);
            break;
        }
        sys->gz(y_new, z_guess, t + step, J, sys->ctx);
        for (int i = 0; i < nz; ++i) {
            delta[i] = -gvec[i];
        }
        if (lin_solve(J, delta, nz)) {
            free(delta);
            free(J);
            free(gvec);
            break;
        }
        double delta_norm = 0.0;
        for (int i = 0; i < nz; ++i) {
            z_guess[i] += delta[i];
            delta_norm += delta[i] * delta[i];
        }
        free(delta);
        free(J);
        free(gvec);
        if (delta_norm < 1e-16) {
            for (int i = 0; i < ny; ++i) {
                y[i] = y_new[i];
            }
            for (int i = 0; i < nz; ++i) {
                z[i] = z_guess[i];
            }
            free(z_guess);
            free(y_new);
            free(dydt);
            return 0;
        }
    }
    free(z_guess);
    free(y_new);
    free(dydt);
    return -1;
}
