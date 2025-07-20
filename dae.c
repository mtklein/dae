#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdlib.h>

_Bool dae_step_euler(struct dae_system const *sys, double step, double *Y, double *Z, double t) {
    int ny = sys->dim_y;
    int nz = sys->dim_z;

    double *dYdt = malloc((size_t)ny * sizeof *dYdt);
    if (!dYdt) {
        return (_Bool)0;
    }
    sys->f(Y, Z, t, dYdt, sys->ctx);

    double *Y_new = malloc((size_t)ny * sizeof *Y_new);
    if (!Y_new) {
        free(dYdt);
        return (_Bool)0;
    }
    for (int i = 0; i < ny; ++i) {
        Y_new[i] = Y[i] + step * dYdt[i];
    }

    double *Z_guess = malloc((size_t)nz * sizeof *Z_guess);
    if (!Z_guess) {
        free(Y_new);
        free(dYdt);
        return (_Bool)0;
    }
    for (int i = 0; i < nz; ++i) {
        Z_guess[i] = Z[i];
    }

    for (int iter = 0; iter < 10; ++iter) {
        double *Gvec = malloc((size_t)nz * sizeof *Gvec);
        if (!Gvec) {
            break;
        }
        sys->g(Y_new, Z_guess, t + step, Gvec, sys->ctx);
        double norm = 0.0;
        for (int i = 0; i < nz; ++i) {
            norm += Gvec[i] * Gvec[i];
        }
        if (norm < 1e-16) {
            free(Gvec);
            for (int i = 0; i < ny; ++i) {
                Y[i] = Y_new[i];
            }
            for (int i = 0; i < nz; ++i) {
                Z[i] = Z_guess[i];
            }
            free(Z_guess);
            free(Y_new);
            free(dYdt);
            return (_Bool)1;
        }
        size_t jm = (size_t)nz * (size_t)nz;
        double *J = malloc(jm * sizeof *J);
        double *Delta = malloc((size_t)nz * sizeof *Delta);
        if (!J || !Delta) {
            free(Delta);
            free(J);
            free(Gvec);
            break;
        }
        sys->gz(Y_new, Z_guess, t + step, J, sys->ctx);
        for (int i = 0; i < nz; ++i) {
            Delta[i] = -Gvec[i];
        }
        if (!lin_solve(J, Delta, nz)) {
            free(Delta);
            free(J);
            free(Gvec);
            break;
        }
        double delta_norm = 0.0;
        for (int i = 0; i < nz; ++i) {
            Z_guess[i] += Delta[i];
            delta_norm += Delta[i] * Delta[i];
        }
        free(Delta);
        free(J);
        free(Gvec);
        if (delta_norm < 1e-16) {
            for (int i = 0; i < ny; ++i) {
                Y[i] = Y_new[i];
            }
            for (int i = 0; i < nz; ++i) {
                Z[i] = Z_guess[i];
            }
            free(Z_guess);
            free(Y_new);
            free(dYdt);
            return (_Bool)1;
        }
    }
    free(Z_guess);
    free(Y_new);
    free(dYdt);
    return (_Bool)0;
}
