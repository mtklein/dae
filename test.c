#include "dae.h"
#include <math.h>
#include <stdio.h>

#define expect(x) \
    if (!(x)) fprintf(stderr, "%s:%d expect(%s)\n", __FILE__,__LINE__,#x), __builtin_debugtrap()

static void expect_close(double got, double expect) {
    expect(fabs(got - expect) <= 1e-3);
}

// First-order decay A -> B at rate k while the total amount A+B stays m.
// The analytic solution A(t)=exp(-k t) lets us check B(t)=m-A(t) after one unit.
struct reaction1 {
    double k;
    double m;
};

static void reaction1_f(double const *y, double const *z, double t, double *dydt, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    (void)z;
    dydt[0] = -r->k * y[0];
}

static void reaction1_g(double const *y, double const *z, double t, double *out, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    out[0] = y[0] + z[0] - r->m;
}

static void reaction1_gz(double const *y, double const *z, double t, double *out, void *ctx) {
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
        expect(dae_step_euler(&sys, 0.001, y, z, t));
        t += 0.001;
    }
    expect_close(y[0], exp(-1.0));
    expect_close(z[0], 1.0 - exp(-1.0));
}

// Second-order reaction A + B -> C with rate k and equal initial A and B.
// With A=B the solution is A(t)=B(t)=1/(1+k t), giving C(t)=2-2*A(t) for comparison.
struct reaction2 {
    double k;
};

static void reaction2_f(double const *y, double const *z, double t, double *dydt, void *ctx) {
    struct reaction2 const *r = ctx;
    (void)t;
    (void)z;
    dydt[0] = -r->k * y[0] * y[1];
    dydt[1] = -r->k * y[0] * y[1];
}

static void reaction2_g(double const *y, double const *z, double t, double *out, void *ctx) {
    (void)t;
    (void)ctx;
    out[0] = y[0] + y[1] + z[0] - 2.0;
}

static void reaction2_gz(double const *y, double const *z, double t, double *out, void *ctx) {
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
        expect(dae_step_euler(&sys, 0.001, y, z, t));
        t += 0.001;
    }
    double expect_a = 1.0 / (1.0 + t);
    expect_close(y[0], expect_a);
    expect_close(y[1], expect_a);
    expect_close(z[0], 2.0 - 2.0 * expect_a);
}

// A pendulum with rigid rod length len and gravitational acceleration g.
// The coordinates (x,y) describe the bob position and (u,v) its velocities.
// A Lagrange multiplier L enforces the length constraint.  Differentiating the
// constraint twice leads to
//   x' = u
//   y' = v
//   u' = L*x
//   v' = L*y - g
//   0  = L*(x^2 + y^2) - (y*g - u*u - v*v)
// The last equation lets the solver compute L so that x^2+y^2 stays len^2.  We
// start slightly offset from vertical and verify the length remains constant.
struct pendulum {
    double len;
    double g;
};

static void pendulum_f(double const *y, double const *z, double t, double *dydt, void *ctx) {
    struct pendulum const *p = ctx;
    (void)t;
    double const L  = z[0];
    double const x  = y[0];
    double const y1 = y[1];
    double const u  = y[2];
    double const v  = y[3];
    dydt[0] = u;
    dydt[1] = v;
    dydt[2] = L * x;
    dydt[3] = L * y1 - p->g;
}

static void pendulum_g(double const *y, double const *z, double t, double *out, void *ctx) {
    struct pendulum const *p = ctx;
    (void)t;
    double const L  = z[0];
    double const x  = y[0];
    double const y1 = y[1];
    double const u  = y[2];
    double const v  = y[3];
    out[0] = L * (x*x + y1*y1) - (y1 * p->g - u*u - v*v);
}

static void pendulum_gz(double const *y, double const *z, double t, double *out, void *ctx) {
    (void)z;
    (void)t;
    (void)ctx;
    double const x  = y[0];
    double const y1 = y[1];
    out[0] = x*x + y1*y1;
}

static void test_pendulum(void) {
    struct pendulum p = {1.0, 1.0};
    struct dae_system sys = {
        4, 1, pendulum_f, pendulum_g, pendulum_gz, &p
    };
    double y[4] = {0.1, sqrt(p.len*p.len - 0.1*0.1), 0.0, 0.0};
    double z[1] = {0.0};
    double t = 0.0;
    for (int i = 0; i < 2000; ++i) {
        expect(dae_step_euler(&sys, 0.0005, y, z, t));
        t += 0.0005;
    }
    expect_close(y[0]*y[0] + y[1]*y[1], p.len * p.len);
}

int main(void) {
    test_reaction1();
    test_reaction2();
    test_pendulum();
    return 0;
}
