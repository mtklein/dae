#ifndef DAE_H
#define DAE_H

struct dae_system {
    int dim_y;
    int dim_z;
    void (*f)(double const *y, double const *z, double t,
              double *dydt, void *ctx);
    void (*g)(double const *y, double const *z, double t,
              double *out, void *ctx);
    void (*gz)(double const *y, double const *z, double t,
               double *out, void *ctx);
    void *ctx;
};

int dae_step_euler(struct dae_system const *sys, double step,
                   double *y, double *z, double t);

#endif
