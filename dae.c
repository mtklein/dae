#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdlib.h>

static _Bool step_euler(struct dae_system const *sys, double t, double dt, double *Y, double *Z,
                        double *dYdt, double *Y_new,
                        double *Z_guess, double *Gvec, double *Delta,
                        double *J) {
    int const nY = sys->dim_y,
              nZ = sys->dim_z;

    sys->f(Y, Z, t, dYdt, sys->ctx);
    for (int i = 0; i < nY; ++i) {
        Y_new[i] = Y[i] + dYdt[i]*dt;
    }
    for (int i = 0; i < nZ; ++i) {
        Z_guess[i] = Z[i];
    }

    for (int iter = 0; iter < 10; ++iter) {
        sys->g(Y_new, Z_guess, t+dt, Gvec, sys->ctx);
        double norm = 0.0;
        for (int i = 0; i < nZ; ++i) {
            norm += Gvec[i] * Gvec[i];
        }
        if (norm < 1e-16) {
            for (int i = 0; i < nY; ++i) {
                Y[i] = Y_new[i];
            }
            for (int i = 0; i < nZ; ++i) {
                Z[i] = Z_guess[i];
            }
            return (_Bool)1;
        }
        sys->gz(Y_new, Z_guess, t+dt, J, sys->ctx);
        for (int i = 0; i < nZ; ++i) {
            Delta[i] = -Gvec[i];
        }
        if (!lin_solve(J, Delta, nZ)) {
            break;
        }
        double delta_norm = 0.0;
        for (int i = 0; i < nZ; ++i) {
            Z_guess[i] += Delta[i];
            delta_norm += Delta[i] * Delta[i];
        }
        if (delta_norm < 1e-16) {
            for (int i = 0; i < nY; ++i) {
                Y[i] = Y_new[i];
            }
            for (int i = 0; i < nZ; ++i) {
                Z[i] = Z_guess[i];
            }
            return (_Bool)1;
        }
    }
    return (_Bool)0;
}

_Bool dae_step_euler(struct dae_system const *sys, double t, double dt, double *Y, double *Z) {
    int const nY = sys->dim_y,
              nZ = sys->dim_z,
              nJ = nZ*nZ;

    double *dYdt    = malloc((size_t)nY * sizeof(double)),
           *Y_new   = malloc((size_t)nY * sizeof(double));
    double *Z_guess = malloc((size_t)nZ * sizeof(double)),
           *Gvec    = malloc((size_t)nZ * sizeof(double)),
           *Delta   = malloc((size_t)nZ * sizeof(double));
    double *J       = malloc((size_t)nJ * sizeof(double));

    _Bool const ok = step_euler(sys, t, dt, Y, Z, dYdt, Y_new, Z_guess, Gvec, Delta, J);

    free(dYdt);
    free(Y_new);
    free(Z_guess);
    free(Gvec);
    free(Delta);
    free(J);
    return ok;
}
