#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdlib.h>

#define MAX_ROUNDS 10
#define CONVERGED  1e-16

#define each(dim) for (int i = 0; i < n##dim; i++)

static _Bool step_euler(struct dae_system const *sys, double t, double dt, double *Y, double *Z,
                        double *dYdt, double *Y_new,
                        double *Z_guess, double *G, double *Delta,
                        double *J) {
    int const nY = sys->dim_y,
              nZ = sys->dim_z;

    sys->f(Y, Z, t, dYdt, sys->ctx);
    each(Y) { Y_new  [i] = Y[i] + dt*dYdt[i]; }
    each(Z) { Z_guess[i] = Z[i]; }

    for (int round = 0; round < MAX_ROUNDS; round++) {
        sys->g(Y_new, Z_guess, t+dt, G, sys->ctx);

        {
            double norm = 0.0;
            each(Z) { norm += G[i] * G[i]; }
            if (norm < CONVERGED) {
                goto success;
            }
        }

        sys->gz(Y_new, Z_guess, t+dt, J, sys->ctx);
        each(Z) { Delta[i] = -G[i]; }
        if (!lin_solve(J, Delta, nZ)) {
            break;
        }

        {
            double norm = 0.0;
            each(Z) {
                Z_guess[i] += Delta[i];
                norm       += Delta[i] * Delta[i];
            }
            if (norm < CONVERGED) {
                goto success;
            }
        }
    }
    return (_Bool)0;

success:
    each(Y) { Y[i] = Y_new[i]; }
    each(Z) { Z[i] = Z_guess[i]; }
    return (_Bool)1;
}

_Bool dae_step_euler(struct dae_system const *sys, double t, double dt, double *Y, double *Z) {
    int const nY = sys->dim_y,
              nZ = sys->dim_z,
              nJ = nZ*nZ;

    double *dYdt    = malloc((size_t)nY * sizeof(double)),
           *Y_new   = malloc((size_t)nY * sizeof(double));
    double *Z_guess = malloc((size_t)nZ * sizeof(double)),
           *G       = malloc((size_t)nZ * sizeof(double)),
           *Delta   = malloc((size_t)nZ * sizeof(double));
    double *J       = malloc((size_t)nJ * sizeof(double));

    _Bool const ok = step_euler(sys, t, dt, Y, Z, dYdt, Y_new, Z_guess, G, Delta, J);

    free(dYdt);
    free(Y_new);
    free(Z_guess);
    free(G);
    free(Delta);
    free(J);
    return ok;
}
