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
    void (*f )(double const *y, double const *z, double t, double *dydt, void *ctx);
    void (*g )(double const *y, double const *z, double t, double  *out, void *ctx);
    void (*gz)(double const *y, double const *z, double t, double  *out, void *ctx);
    void *ctx;
};

// Advance the system by one Euler step and update z using Newton iteration.
_Bool dae_step_euler(struct dae_system const *sys, double step, double *y, double *z, double t);
