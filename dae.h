#pragma once

// Represents a differential-algebraic system
//     dY/dt = f(Y, Z, t)
//         0 = g(Y, Z, t)
// dim_Y and dim_Z give the counts of Y and Z values.
// f evaluates the differential equation, g the algebraic constraint,
// and gz is the derivative of g with respect to Z.
// ctx is user data forwarded to the callbacks.
struct dae_system {
    int dim_y,
        dim_z;
    void (*f )(double const *Y, double const *Z, double t, double *dYdt, void *ctx);
    void (*g )(double const *Y, double const *Z, double t, double *G   , void *ctx);
    void (*gz)(double const *Y, double const *Z, double t, double *dGdZ, void *ctx);
    void *ctx;
};

// Advance the system by one Euler step and update z using Newton iteration.
_Bool dae_step_euler(struct dae_system const *sys, double step, double *Y, double *Z, double t);
