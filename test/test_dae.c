#include "dae.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void expect_close(double got, double expect) {
    if (fabs(got - expect) > 1e-3) {
        fprintf(stderr, "expect %f got %f\n", expect, got);
        abort();
    }
}

struct reaction1 {
    double k;
    double m;
};

static void reaction1_f(double const *y, double const *z, double t,
                        double *dydt, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    (void)z;
    dydt[0] = -r->k * y[0];
}

static void reaction1_g(double const *y, double const *z, double t,
                        double *out, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    out[0] = y[0] + z[0] - r->m;
}

static void reaction1_gz(double const *y, double const *z, double t,
                         double *out, void *ctx) {
    (void)y;
    (void)z;
    (void)t;
    (void)ctx;
    out[0] = 1.0;
}

static void test_reaction1(void) {
    struct reaction1 r = {1.0, 1.0};
    struct dae_system sys = {
        1, 1, reaction1_f, reaction1_g, reaction1_gz, &r
    };
    double y[1] = {1.0};
    double z[1] = {0.0};
    double t = 0.0;
    for (int i = 0; i < 1000; ++i) {
        if (dae_step_euler(&sys, 0.001, y, z, t)) {
            abort();
        }
        t += 0.001;
    }
    expect_close(y[0], exp(-1.0));
    expect_close(z[0], 1.0 - exp(-1.0));
}

struct reaction2 {
    double k;
};

static void reaction2_f(double const *y, double const *z, double t,
                        double *dydt, void *ctx) {
    struct reaction2 const *r = ctx;
    (void)t;
    (void)z;
    dydt[0] = -r->k * y[0] * y[1];
    dydt[1] = -r->k * y[0] * y[1];
}

static void reaction2_g(double const *y, double const *z, double t,
                        double *out, void *ctx) {
    (void)t;
    (void)ctx;
    out[0] = y[0] + y[1] + z[0] - 2.0;
}

static void reaction2_gz(double const *y, double const *z, double t,
                         double *out, void *ctx) {
    (void)y;
    (void)z;
    (void)t;
    (void)ctx;
    out[0] = 1.0;
}

static void test_reaction2(void) {
    struct reaction2 r = {1.0};
    struct dae_system sys = {
        2, 1, reaction2_f, reaction2_g, reaction2_gz, &r
    };
    double y[2] = {1.0, 1.0};
    double z[1] = {0.0};
    double t = 0.0;
    for (int i = 0; i < 1000; ++i) {
        if (dae_step_euler(&sys, 0.001, y, z, t)) {
            abort();
        }
        t += 0.001;
    }
    double expect_a = 1.0 / (1.0 + t);
    expect_close(y[0], expect_a);
    expect_close(y[1], expect_a);
    expect_close(z[0], 2.0 - 2.0 * expect_a);
}

int main(void) {
    test_reaction1();
    test_reaction2();
    return 0;
}
