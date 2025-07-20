#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdio.h>

#define len(x) (int)(sizeof x / sizeof *x)

#define expect(x) \
    if (!(x)) fprintf(stderr, "%s:%d expect(%s)\n", __FILE__,__LINE__,#x), __builtin_debugtrap()

static void expect_equiv(double x, double y) {
    expect((x <= y && y <= x) ||
           (x != x && y != y));
}
static void expect_close(double x, double y) {
    expect(fabs(x - y) <= 1e-3);
}

static void test_lin_solve_2x2(void) {
    // Solve the row-major system
    //   2x + 3y = 5
    //    x +  y = 2
    double M[] = {
        2.0, 3.0,
        1.0, 1.0,
    };
    double V[] = {5.0, 2.0};
    expect(lin_solve(M, V, len(V)));
    expect_equiv(V[0], 1.0);
    expect_equiv(V[1], 1.0);
}

static void test_lin_solve_pivot(void) {
    // Solve the row-major system requiring a pivot swap
    //        2y + 3z = 1
    //    x + 2y + 3z = 2
    //   4x + 5y + 6z = 3
    double M[] = {
        0.0, 2.0, 3.0,
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
    };
    double V[] = {1.0, 2.0, 3.0};
    expect(lin_solve(M, V, len(V)));
    expect_equiv(V[0],    1.0);
    expect_equiv(V[1],   -3.0);
    expect_equiv(V[2],  7/3.0);
}

// First-order decay A -> B at rate k while the total amount A+B stays m.
// The analytic solution A(t)=exp(-k t) lets us check B(t)=m-A(t) after one unit.
struct reaction1 {
    double const k,m;
};

static void reaction1_f(double const *Y, double const *Z, double t, double *dYdt, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    (void)Z;
    double const A = Y[0];
    dYdt[0] = -r->k * A;
}

static void reaction1_g(double const *Y, double const *Z, double t, double *G, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    double const A = Y[0],
                 B = Z[0];
    G[0] = A + B - r->m;
}

static void reaction1_gz(double const *Y, double const *Z, double t, double *dGdZ, void *ctx) {
    (void)Y;
    (void)Z;
    (void)t;
    (void)ctx;
    dGdZ[0] = 1.0;
}

static void test_reaction1(void) {
    struct reaction1 r = {.k=1.0, .m=1.0};
    double Y[] = {1.0};
    double Z[] = {0.0};
    double t   = 0.0;
    struct dae_system sys = {len(Y), len(Z), reaction1_f, reaction1_g, reaction1_gz, &r};
    for (int i = 0; i < 1000; ++i) {
        expect(dae_step_euler(&sys, 0.001, Y, Z, t));
        t += 0.001;
    }
    expect_close(Y[0], exp(-1.0));
    expect_equiv(Z[0], 1.0 - Y[0]);
}

// Second-order reaction A + B -> C with rate k and equal initial A and B.
// With A=B the solution is A(t)=B(t)=1/(1+k t), giving C(t)=2-2*A(t) for comparison.
struct reaction2 {
    double const k;
};

static void reaction2_f(double const *Y, double const *Z, double t, double *dYdt, void *ctx) {
    struct reaction2 const *r = ctx;
    (void)t;
    (void)Z;
    double const A = Y[0],
                 B = Y[1];
    dYdt[0] = -r->k * A * B;
    dYdt[1] = -r->k * A * B;
}

static void reaction2_g(double const *Y, double const *Z, double t, double *G, void *ctx) {
    (void)t;
    (void)ctx;
    double const A = Y[0],
                 B = Y[1],
                 C = Z[0];
    G[0] = A + B + C - 2.0;
}

static void reaction2_gz(double const *Y, double const *Z, double t, double *dGdZ, void *ctx) {
    (void)Y;
    (void)Z;
    (void)t;
    (void)ctx;
    dGdZ[0] = 1.0;
}

static void test_reaction2(void) {
    struct reaction2 r = {.k=1.0};
    double Y[] = {1.0, 1.0};
    double Z[] = {0.0};
    double t   = 0.0;
    struct dae_system sys = {len(Y), len(Z), reaction2_f, reaction2_g, reaction2_gz, &r};

    for (int i = 0; i < 1000; ++i) {
        expect(dae_step_euler(&sys, 0.001, Y, Z, t));
        t += 0.001;
    }
    double const expect_a = 1.0 / (1.0 + t);
    expect_equiv(Y[0], Y[1]);
    expect_close(Y[0], expect_a);
    expect_close(Y[1], expect_a);
    expect_close(Z[0], 2.0 - 2.0 * expect_a);
}

// A pendulum with rigid rod length len and gravitational acceleration g.
// The coordinates (x,y) describe the bob position and (u,v) its velocities.
// A Lagrange multiplier l enforces the length constraint.  Differentiating the
// constraint twice leads to
//   x' = u
//   y' = v
//   u' = l*x
//   v' = l*y - g
//   0  = l*(x^2 + y^2) - (y*g - u*u - v*v)
// The last equation lets the solver compute l so that x^2+y^2 stays len^2.  We
// start slightly offset from vertical and verify the length remains constant.
struct pendulum {
    double len;
    double g;
};

static void pendulum_f(double const *Y, double const *Z, double t, double *dYdt, void *ctx) {
    struct pendulum const *p = ctx;
    (void)t;
    double const x = Y[0],
                 y = Y[1],
                 u = Y[2],
                 v = Y[3],
                 l = Z[0];
    dYdt[0] = u;
    dYdt[1] = v;
    dYdt[2] = l * x;
    dYdt[3] = l * y - p->g;
}

static void pendulum_g(double const *Y, double const *Z, double t, double *G, void *ctx) {
    struct pendulum const *p = ctx;
    (void)t;
    double const x = Y[0],
                 y = Y[1],
                 u = Y[2],
                 v = Y[3],
                 l = Z[0];
    G[0] = l * (x*x + y*y) - (y * p->g - u*u - v*v);
}

static void pendulum_gz(double const *Y, double const *Z, double t, double *dGdZ, void *ctx) {
    (void)Z;
    (void)t;
    (void)ctx;
    double const x = Y[0],
                 y = Y[1];
    dGdZ[0] = x*x + y*y;
}

static void test_pendulum(void) {
    struct pendulum p = {1.0, 1.0};
    double Y[] = {0.1, sqrt(p.len*p.len - 0.1*0.1), 0.0, 0.0};
    double Z[] = {0.0};
    double t   = 0.0;
    struct dae_system sys = {len(Y), len(Z), pendulum_f, pendulum_g, pendulum_gz, &p};

    for (int i = 0; i < 2000; ++i) {
        expect(dae_step_euler(&sys, 0.0005, Y, Z, t));
        t += 0.0005;
    }
    expect_close(Y[0]*Y[0] + Y[1]*Y[1], p.len*p.len);
}

int main(void) {
    test_lin_solve_2x2();
    test_lin_solve_pivot();
    test_reaction1();
    test_reaction2();
    test_pendulum();
    return 0;
}
