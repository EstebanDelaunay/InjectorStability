/* =======================================================================================================

CC99='mpicc -std=c99' qcc -Wall -O2 -autolink -D_MPI=1 main.c -o main -lm -lfb_tiny ; mpirun -np 8 ./main

qcc -autolink -Wall -O2  main.c -o main -lm -lfb_tiny ; ./main

======================================================================================================= */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"

#include "view.h"
#include "common.h"
#include "fractions.h"
#include "draw.h"


// =======================================================
// Time parameters =======================================

const double tEnd = 1e0;
const double tStep = 1e-1;


// =======================================================
// Space parameters ======================================

const double jetThickness = 1e0;
const double domainLength = 50. * jetThickness;


// =======================================================
// Physical parameters ===================================
// L for liquid ; G for gaz

scalar f0[];
char name[80];

// Density
const double rhoL = 1., rhoG = 1e-3;
// Viscosity
const double muL= 1e-3, muG = 1e-5;
// Velocity
const double u0 = 1e-1;
// Gravity
const double gravity = 10;
// Surface tension
const double sigma = 71.97e-3;

const double Re = rhoL * u0 * jetThickness / muL;
const double Fr = u0 * u0 / (gravity * jetThickness);
const double We = rhoL * u0 * u0 * jetThickness / sigma;


// =======================================================
// Mesh parameters =======================================

const int maxlevel = 9;
const double uemax = 1e-3;


// =======================================================
// Boundary conditions ===================================

// Inflow

u.n[left] = dirichlet(f0[] * u0);
u.t[left] = dirichlet(0.);
p[left]   = neumann(0.);
f[left]   = f0[];

// Outflow

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

f[right]   = dirichlet(0.);


// =======================================================
// main ==================================================
int main()
{
    periodic(top);

    //init_grid(1 << maxlevel);
    init_grid(1 << 7);

    origin(0., -domainLength/2.);
    size(domainLength);

    rho1 = rhoL;
    rho2 = rhoG;

    mu1 = muL;
    mu2 = muG;

    f.sigma = sigma;

    run();
}


// =======================================================
// Initial conditions ====================================
event init(t = 0)
{
    TOLERANCE = 1e-3 [*];

    G.x = gravity;

    /* static refinement down to *maxlevel* in a plan twice the width. */
    refine(sq(y) < 2. * sq(jetThickness) && level < maxlevel);

    fraction(f0, sq(jetThickness) - sq(y));

    f0.refine = f0.prolongation = fraction_refine;
    restriction({f0}); // for boundary conditions on levels

    foreach ()
    {
        f[] = f0[];
        u.x[] = u0 * f[];
    }
}


// =======================================================
// Log ===================================================
event logfile(i++)
{
    if (i == 0)
    {
        fprintf(stderr, "-----------------------------------------------------\n");
        fprintf(stderr, "   Re   |    Fr   |    We   |    u0   | gravity |  sigma  \n");
        fprintf(stderr, "%1.1e | %1.1e | %1.1e | %1.1e | %1.1e | %1.1e    \n", Re, Fr, We, u0, gravity, sigma);
        fprintf(stderr, "-----------------------------------------------------\n");
        fprintf(stderr, "   t    |    dt    | p.i | f.i | u.i | ncells | perf \n");
        fprintf(stderr, "        |          |     |     |     |        |      \n");
    }

    fprintf(stderr, "%.5f | %.2e |  %2d |  %2d |  %2d | %6ld | %.2f\n", t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t);
}


// =======================================================
// Movie =================================================
event movie(t += tStep; t <= tEnd)
{ 
    scalar omega[];
    vorticity(u, omega);
    view(fov = 25., quat = {0., 0., cos(-pi/4.), cos(pi/4.)}, tx = 0., ty = 0.5);

    clear();
    draw_vof("f");
    squares("f", min=0., max=1.);

    sprintf(name, "out_p_%d_%.3f_%.3f_%.3f__%.5f.dat", LEVEL, We, Ma, A, f.sigma);
    draw_string(str = "Re", pos = 3, size = 25, lc = {255.,255.,255.}, lw = 4);
    save("movie.mp4");
}


// =======================================================
// Mesh adaptation =======================================

event adapt(i++)
{
    adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax, uemax}, maxlevel, 5);
}
