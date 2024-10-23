/**
# Atomisation of a pulsed liquid jet

CC99='mpicc -std=c99' qcc -Wall -O2 -autolink -D_MPI=1 main.c -o main -lm -lfb_tiny ; mpirun -np 8 ./main

*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"

#include "view.h"
#include "common.h"
#include "fractions.h"



// =======================================================
// Time parameters =======================================
const double tEnd = 5e0;
const double tStep = 1e-1;



// =======================================================
// Space parameters ======================================
const double jetThickness = 0.1;
const double domainLength = 100. * jetThickness;



// =======================================================
// Physical parameters ===================================
scalar f0[];
const double Re = 5;
const double Fr = 1.;

const double SIGMA = 3e-5;
const double u0 = 0.1;



// =======================================================
// Mesh parameters =======================================
const int maxlevel = 8; // default maximum level of refinement
const double uemax = 0.1; // error threshold on velocity


// The inflow condition fixes the velocity
u.n[left] = dirichlet(f0[] * u0);
u.t[left] = dirichlet(0);
p[left] = neumann(0);

// Outflow uses standard Neumann/Dirichlet conditions.
u.n[right] = neumann(0);
u.t[right] = neumann(0);
p[right]    = dirichlet(0.);
pf[right]   = dirichlet(0.);

// Fraction boundaries
f[left] = f0[];
f[right] = dirichlet(0);



// =======================================================
// main ==================================================
int main()
{
    //periodic(front);

    init_grid(1 << maxlevel);
    origin(0., -domainLength/2.);
    size(domainLength);

    rho1 = 1. [0];
    rho2 = rho1 / 815.;

    mu1 = 2. * u0 * jetThickness / Re * rho1;
    mu2 = 2. * u0 * jetThickness / Re * rho2;

    f.sigma = SIGMA;

    run();
}



// =======================================================
// Initial conditions ====================================
event init(t = 0)
{
    TOLERANCE = 1e-4 [*];

    G.x = sq(u0 / Fr) / jetThickness; // maybe domainLength

    /* We use a static refinement down to *maxlevel* in a plan twice the width. */
    refine(fabs(y) < 2. * jetThickness && level < maxlevel);

    fraction(f0, fabs(y) < jetThickness );

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
        fprintf(stderr, "t dt mgp.i mgpf.i mgu.i grid->tn perf.t \n");
    }

    fprintf(stderr, "%.5f %.2e %2d %2d %2d %ld %.2f\n",
            t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t);
}



// =======================================================
// Movie =================================================

event movie(t += tStep; t <= tEnd)
{ 
    scalar omega[];
    vorticity(u, omega);
    view(tx = -.5);
    clear();
    draw_vof("f");
    //squares("omega", linear = true, spread = 10);
    squares("f", min=0., max=1.);
    //cells();
    box();
    save("movie.mp4");
}



// =======================================================
// Mesh adaptation =======================================

event adapt(i++)
{
    adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax, uemax}, maxlevel, 3);
}
