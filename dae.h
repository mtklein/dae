#pragma once

// TODO: replace this comment with a short paragraph describing what dae_system represents,
//       and a brief description of each of its parts if not obvious from that description.
struct dae_system {
    int dim_y,
        dim_z;
    void (*f )(double const *y, double const *z, double t, double *dydt, void *ctx);
    void (*g )(double const *y, double const *z, double t, double  *out, void *ctx);
    void (*gz)(double const *y, double const *z, double t, double  *out, void *ctx);
    void *ctx;
};

// TODO: replace this comment with a one-line description of what this function does.
_Bool dae_step_euler(struct dae_system const *sys, double step, double *y, double *z, double t);
