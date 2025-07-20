#include "dae.h"
#include "lin.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

static _Bool plot_enabled;

static void plot_comment(char const *name) {
    if (plot_enabled) {
        printf("# %s\n", name);
    }
}

static void plot_line(double t, double const *Y, int nY,
                      double const *Z, int nZ) {
    if (!plot_enabled) {
        return;
    }
    printf("%g", t);
    for (int i = 0; i < nY; ++i) {
        printf(" %g", Y[i]);
    }
    for (int i = 0; i < nZ; ++i) {
        printf(" %g", Z[i]);
    }
    putchar('\n');
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

static void reaction1_f(double *dYdt, double t, double const *Y, double const *Z, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    (void)Z;
    double const A = Y[0];
    dYdt[0] = -r->k * A;
}

static void reaction1_g(double *G, double t, double const *Y, double const *Z, void *ctx) {
    struct reaction1 const *r = ctx;
    (void)t;
    double const A = Y[0],
                 B = Z[0];
    G[0] = A + B - r->m;
}

static void reaction1_gz(double *dGdZ, double t, double const *Y, double const *Z, void *ctx) {
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
    plot_comment("reaction1");
    for (int i = 0; i < 1000; ++i) {
        double const dt = 1.0/1000;
        expect(dae_step_euler(&sys, t, dt, Y, Z));
        t += dt;
        plot_line(t, Y, len(Y), Z, len(Z));
    }
    expect_close(Y[0], exp(-1.0));
    expect_equiv(Z[0], 1.0 - Y[0]);
}

// Second-order reaction A + B -> C with rate k and equal initial A and B.
// With A=B the solution is A(t)=B(t)=1/(1+k t), giving C(t)=2-2*A(t) for comparison.
struct reaction2 {
    double const k;
};

static void reaction2_f(double *dYdt, double t, double const *Y, double const *Z, void *ctx) {
    struct reaction2 const *r = ctx;
    (void)t;
    (void)Z;
    double const A = Y[0],
                 B = Y[1];
    dYdt[0] = -r->k * A * B;
    dYdt[1] = -r->k * A * B;
}

static void reaction2_g(double *G, double t, double const *Y, double const *Z, void *ctx) {
    (void)t;
    (void)ctx;
    double const A = Y[0],
                 B = Y[1],
                 C = Z[0];
    G[0] = A + B + C - 2.0;
}

static void reaction2_gz(double *dGdZ, double t, double const *Y, double const *Z, void *ctx) {
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

    plot_comment("reaction2");

    for (int i = 0; i < 1000; ++i) {
        double const dt = 1.0/1000;
        expect(dae_step_euler(&sys, t, dt, Y, Z));
        t += dt;
        plot_line(t, Y, len(Y), Z, len(Z));
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

static void pendulum_f(double *dYdt, double t, double const *Y, double const *Z, void *ctx) {
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

static void pendulum_g(double *G, double t, double const *Y, double const *Z, void *ctx) {
    struct pendulum const *p = ctx;
    (void)t;
    double const x = Y[0],
                 y = Y[1],
                 u = Y[2],
                 v = Y[3],
                 l = Z[0];
    G[0] = l * (x*x + y*y) - (y * p->g - u*u - v*v);
}

static void pendulum_gz(double *dGdZ, double t, double const *Y, double const *Z, void *ctx) {
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

    plot_comment("pendulum");

    for (int i = 0; i < 2000; ++i) {
        double const dt = 1.0/2000;
        expect(dae_step_euler(&sys, t, dt, Y, Z));
        t += dt;
        plot_line(t, Y, len(Y), Z, len(Z));
    }
    expect_close(Y[0]*Y[0] + Y[1]*Y[1], p.len*p.len);
}

// Two independent pendulums sharing the same gravity but with different rod
// lengths.  Each pendulum uses its own Lagrange multiplier to enforce the
// length constraint.
struct pendulum_pair {
    double len1,len2,g;
};

static void pendulum_pair_f(double *dYdt, double t, double const *Y, double const *Z, void *ctx) {
    struct pendulum_pair const *p = ctx;
    (void)t;
    double const x1 = Y[0], y1 = Y[1], u1 = Y[2], v1 = Y[3];
    double const x2 = Y[4], y2 = Y[5], u2 = Y[6], v2 = Y[7];
    double const l1 = Z[0], l2 = Z[1];
    dYdt[0] = u1;
    dYdt[1] = v1;
    dYdt[2] = l1 * x1;
    dYdt[3] = l1 * y1 - p->g;
    dYdt[4] = u2;
    dYdt[5] = v2;
    dYdt[6] = l2 * x2;
    dYdt[7] = l2 * y2 - p->g;
}

static void pendulum_pair_g(double *G, double t, double const *Y, double const *Z, void *ctx) {
    struct pendulum_pair const *p = ctx;
    (void)t;
    double const x1 = Y[0], y1 = Y[1], u1 = Y[2], v1 = Y[3], l1 = Z[0];
    double const x2 = Y[4], y2 = Y[5], u2 = Y[6], v2 = Y[7], l2 = Z[1];
    G[0] = l1 * (x1*x1 + y1*y1) - (y1 * p->g - u1*u1 - v1*v1);
    G[1] = l2 * (x2*x2 + y2*y2) - (y2 * p->g - u2*u2 - v2*v2);
}

static void pendulum_pair_gz(double *dGdZ, double t, double const *Y, double const *Z, void *ctx) {
    (void)t;
    (void)Z;
    (void)ctx;
    double const x1 = Y[0], y1 = Y[1];
    double const x2 = Y[4], y2 = Y[5];
    dGdZ[0*2 + 0] = x1*x1 + y1*y1;
    dGdZ[0*2 + 1] = 0.0;
    dGdZ[1*2 + 0] = 0.0;
    dGdZ[1*2 + 1] = x2*x2 + y2*y2;
}

static void test_pendulum_pair(void) {
    struct pendulum_pair p = {1.0, 1.5, 1.0};
    double Y[] = {
        0.1, sqrt(p.len1*p.len1 - 0.1*0.1), 0.0, 0.0,
        0.2, sqrt(p.len2*p.len2 - 0.2*0.2), 0.0, 0.0,
    };
    double Z[] = {0.0, 0.0};
    double t   = 0.0;
    struct dae_system sys = {len(Y),len(Z), pendulum_pair_f,pendulum_pair_g,pendulum_pair_gz, &p};

    plot_comment("pendulum_pair");

    for (int i = 0; i < 4000; ++i) {
        double const dt = 1.0/4000;
        expect(dae_step_euler(&sys, t, dt, Y, Z));
        t += dt;
        plot_line(t, Y, len(Y), Z, len(Z));
    }
    expect_close(Y[0]*Y[0] + Y[1]*Y[1], p.len1*p.len1);
    expect_close(Y[4]*Y[4] + Y[5]*Y[5], p.len2*p.len2);
}

int main(void) {
    plot_enabled = getenv("PLOT") != NULL;
    test_lin_solve_2x2();
    test_lin_solve_pivot();
    test_reaction1();
    test_reaction2();
    test_pendulum();
    test_pendulum_pair();
    return 0;
}
