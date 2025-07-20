#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdlib.h>

#define MAX_ROUNDS 10
#define CONVERGED  (1e-8 * 1e-8)

#define each(dim) for (int i = 0; i < n##dim; i++)

static _Bool step_euler(struct dae_system const *sys, double t, double dt, double *Y, double *Z,
                        double *Y_scratch, double *Z_guess, double *G, double *J) {
    int const nY = sys->dim_y,
              nZ = sys->dim_z;

    double *dYdt  = Y_scratch,
           *Y_new = Y_scratch;

    sys->f(dYdt, t,Y,Z,sys->ctx);
    each(Y) { Y_new  [i] = Y[i] + dt*dYdt[i]; }
    each(Z) { Z_guess[i] = Z[i]; }

    for (int rounds = MAX_ROUNDS; rounds --> 0;) {
        sys->g(G, t+dt,Y_new,Z_guess,sys->ctx);

        {
            double norm = 0.0;
            each(Z) {
                norm += G[i] * G[i];
            }
            if (norm < CONVERGED) {
                goto success;
            }
        }

        sys->gz(J, t+dt,Y_new,Z_guess,sys->ctx);
        each(Z) { G[i] = -G[i]; }
        if (!lin_solve(J, G, nZ)) {
            break;
        }

        {
            double norm = 0.0;
            each(Z) {
                norm       += G[i] * G[i];
                Z_guess[i] += G[i];
            }
            if (norm < CONVERGED) {
                goto success;
            }
        }
    }
    return (_Bool)0;

success:
    each(Y) { Y[i] = Y_new  [i]; }
    each(Z) { Z[i] = Z_guess[i]; }
    return (_Bool)1;
}

_Bool dae_step_euler(struct dae_system const *sys, double t, double dt, double *Y, double *Z) {
    int const nY = sys->dim_y,
              nZ = sys->dim_z,
              nJ = nZ*nZ;

    double *Y_scratch = malloc((size_t)nY * sizeof(double));
    double *Z_scratch = malloc((size_t)nZ * sizeof(double)),
           *G         = malloc((size_t)nZ * sizeof(double));
    double *J         = malloc((size_t)nJ * sizeof(double));

    _Bool const ok = step_euler(sys, t, dt, Y, Z, Y_scratch, Z_scratch, G, J);

    free(Y_scratch);
    free(Z_scratch);
    free(G);
    free(J);
    return ok;
}
