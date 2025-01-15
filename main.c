/* =======================================================================================================

CC99='mpicc -std=c99' qcc -Wall -O2 -autolink -D_MPI=1 main.c -o main -lm -lfb_tiny ; mpirun -np 8 ./main

qcc -autolink -Wall -O2  main.c -o main -lm -lfb_tiny ; ./main

======================================================================================================= */

#include "grid/multigrid.h"

#include "navier-stokes/centered.h"

#include "two-phase.h"
#include "tension.h"
#include "reduced.h"

#include "navier-stokes/conserving.h"

#include "view.h"
#include "common.h"
#include "draw.h"

//#include "navier-stokes/perfs.h"

// =======================================================
// Time parameters =======================================

const double tEnd = 5e-1;
const double tStep = 2e-3;


// =======================================================
// Space parameters ======================================

const double jetThickness = 3e-3; // 5 or 10mm
const double domainLength = 50. * jetThickness;


// =======================================================
// Physical parameters ===================================
// L for liquid ; G for gaz

scalar f0[], U0[];
char name[80];

// Density
const double rhoL = 1000., rhoG = 1.225;
// Viscosity
const double muL= 1e-3, muG = 1.8e-5;
// Velocity
const double u0 = 5e-1; 
// Gravity.
const double gravity = 9.81;
// Surface tension
const double sigma = 72e-3;

const double Re = rhoL * u0 * jetThickness / muL; // 1800
const double Fr = u0 * u0 / (gravity * jetThickness);
const double We = rhoL * u0 * u0 * jetThickness / sigma; // fixed Weber number to find speed with jetThickness (We = 0.002 -> 6k, smooth around = ???)
//double lambda = sqrt(sigma / ((rhoL-rhoG)*gravity)) * 1e3; // 2.7mm


// =======================================================
// Mesh parameters =======================================

const int maxlevel = 8;
const double uemax = 1e-3;


// =======================================================
// Boundary conditions ===================================

// Inflow

u.n[left] = dirichlet(f0[] * U0[]);
u.t[left] = dirichlet(0.);
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
    //periodic(front);

    init_grid(1 << maxlevel);

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
    //refine(sq(y) < 2. * sq(jetThickness) && level < maxlevel);

    fraction(f0, sq(jetThickness) - sq(y));

    //f0.refine = f0.prolongation = fraction_refine;
    restriction({f0}); // for boundary conditions on levels

    foreach ()
    {
        //f[] = f0[];
        U0[] = u0 * (1 - sq(y/jetThickness));
        u.x[] = U0[] * f[];
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
        fprintf(stderr, " t / t_max | dt | p.i | f.i | u.i | ncells | perf \n");
        fprintf(stderr, "--------------------------------------------------------\n");

    }

    fprintf(stderr, "%.5f / %.5f | %.2e |  %2d |  %2d |  %2d | %6ld | %.2f\n", t, tEnd, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t);
}


// =======================================================
// Movie =================================================

event movie(t += tStep; t <= tEnd)
{ 
    view(fov = 20., quat = {0., 0., cos(-pi/4.), cos(pi/4.)}, tx = 0., ty = 0.5, width = 750, height = 750);

    // Figure for f
    clear();
    draw_vof("f");
    draw_vof ("f", filled = 1, fc = {0.1, 0.1, 0.8});
    //begin_mirror({0,-1});
    //draw_vof ("f", filled = 1, fc = {0.1, 0.1, 0.8});
    //end_mirror();
    sprintf(name, "movie_f_%1.1e.mp4", jetThickness);
    save(name);

    // Figure for u.x
    clear();
    draw_vof("f");
    squares("u.x", linear = true);
    //begin_mirror({0,-1});
    //draw_vof("f");
    //squares("u.x", linear = true);
    //end_mirror();
    sprintf(name, "movie_u_%1.1e.mp4", jetThickness);
    save(name);

    // Figure for the vorticity
    scalar omega[];
    vorticity(u, omega);
    clear();
    draw_vof("f");
    squares("omega", linear = true);
    //begin_mirror({0,-1});
    //squares("omega", linear = true);
    //end_mirror();
    sprintf(name, "movie_w_%1.1e.mp4", jetThickness);
    save(name);
}


// =======================================================
// Mesh adaptation =======================================

// event adapt(i++)
// {
//     adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax, uemax}, maxlevel, 5);
// }
