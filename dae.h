#pragma once

// Represents a differential-algebraic system
//     dy/dt = f(y, z, t)
//         0 = g(y, z, t)
// dim_y and dim_z give the counts of y and z values.
// f evaluates the differential equation, g the algebraic constraint,
// and gz is the derivative of g with respect to z.
// ctx is user data forwarded to the callbacks.
struct dae_system {
    int dim_y,
        dim_z;
    void (*f )(double const *Y, double const *Z, double t, double *Dydt, void *ctx);
    void (*g )(double const *Y, double const *Z, double t, double *Out, void *ctx);
    void (*gz)(double const *Y, double const *Z, double t, double *Out, void *ctx);
    void *ctx;
};

// Advance the system by one Euler step and update z using Newton iteration.
_Bool dae_step_euler(struct dae_system const *sys, double step, double *Y, double *Z, double t);
