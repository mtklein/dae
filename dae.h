#pragma once

// Represents a differential-algebraic system
//     dY/dt = f(t,Y,Z)
//         0 = g(t,Y,Z)
// dim_Y and dim_Z give the counts of Y and Z values.
// f evaluates the differential equation, g the algebraic constraint,
// and gz is the derivative of g with respect to Z.
// ctx is user data forwarded to the callbacks.
struct dae_system {
    int dim_y,
        dim_z;
    void (*f )(double *dYdt, double t, double const *Y, double const *Z, void *ctx);
    void (*g )(double *G   , double t, double const *Y, double const *Z, void *ctx);
    void (*gz)(double *dGdZ, double t, double const *Y, double const *Z, void *ctx);
    void *ctx;
};

// Advance the system by one Euler step and update z using Newton iteration.
_Bool dae_step_euler(struct dae_system const *sys, double t, double dt, double *Y, double *Z);
