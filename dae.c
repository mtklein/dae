#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdlib.h>

static _Bool dae_step_euler_inner(struct dae_system const *sys, double step,
                                  double *Y, double *Z, double t,
                                  double *dYdt, double *Y_new,
                                  double *Z_guess, double *Gvec,
                                  double *J, double *Delta);

_Bool dae_step_euler(struct dae_system const *sys, double step, double *Y, double *Z, double t) {
    int ny = sys->dim_y;
    int nz = sys->dim_z;

    size_t jm = (size_t)nz * (size_t)nz;
    double *dYdt  = malloc((size_t)ny * sizeof *dYdt);
    double *Y_new = malloc((size_t)ny * sizeof *Y_new);
    double *Z_guess = malloc((size_t)nz * sizeof *Z_guess);
    double *Gvec  = malloc((size_t)nz * sizeof *Gvec);
    double *J     = malloc(jm * sizeof *J);
    double *Delta = malloc((size_t)nz * sizeof *Delta);
    if (!Delta || !J || !Gvec || !Z_guess || !Y_new || !dYdt) {
        free(Delta);
        free(J);
        free(Gvec);
        free(Z_guess);
        free(Y_new);
        free(dYdt);
        return (_Bool)0;
    }

    _Bool const ok = dae_step_euler_inner(sys, step, Y, Z, t, dYdt, Y_new,
                                          Z_guess, Gvec, J, Delta);

    free(Delta);
    free(J);
    free(Gvec);
    free(Z_guess);
    free(Y_new);
    free(dYdt);
    return ok;
}

static _Bool dae_step_euler_inner(struct dae_system const *sys, double step,
                                  double *Y, double *Z, double t,
                                  double *dYdt, double *Y_new,
                                  double *Z_guess, double *Gvec,
                                  double *J, double *Delta) {
    int ny = sys->dim_y;
    int nz = sys->dim_z;

    sys->f(Y, Z, t, dYdt, sys->ctx);
    for (int i = 0; i < ny; ++i) {
        Y_new[i] = Y[i] + step * dYdt[i];
    }
    for (int i = 0; i < nz; ++i) {
        Z_guess[i] = Z[i];
    }

    for (int iter = 0; iter < 10; ++iter) {
        sys->g(Y_new, Z_guess, t + step, Gvec, sys->ctx);
        double norm = 0.0;
        for (int i = 0; i < nz; ++i) {
            norm += Gvec[i] * Gvec[i];
        }
        if (norm < 1e-16) {
            for (int i = 0; i < ny; ++i) {
                Y[i] = Y_new[i];
            }
            for (int i = 0; i < nz; ++i) {
                Z[i] = Z_guess[i];
            }
            return (_Bool)1;
        }
        sys->gz(Y_new, Z_guess, t + step, J, sys->ctx);
        for (int i = 0; i < nz; ++i) {
            Delta[i] = -Gvec[i];
        }
        if (!lin_solve(J, Delta, nz)) {
            break;
        }
        double delta_norm = 0.0;
        for (int i = 0; i < nz; ++i) {
            Z_guess[i] += Delta[i];
            delta_norm += Delta[i] * Delta[i];
        }
        if (delta_norm < 1e-16) {
            for (int i = 0; i < ny; ++i) {
                Y[i] = Y_new[i];
            }
            for (int i = 0; i < nz; ++i) {
                Z[i] = Z_guess[i];
            }
            return (_Bool)1;
        }
    }
    return (_Bool)0;
}
