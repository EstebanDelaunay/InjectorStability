/* =======================================================================================================

CC99='mpicc -std=c99' qcc -Wall -O2 -autolink -D_MPI=1 main.c -o main -lm -lfb_tiny ; mpirun -np 8 ./main

qcc -autolink -Wall -O2  main.c -o main -lm -lfb_tiny ; ./main

======================================================================================================= */

//#include "grid/multigrid.h"

#include "navier-stokes/centered.h"

#include "two-phase.h"
#include "tension.h"
#include "reduced.h"

#include "navier-stokes/conserving.h"

#include "view.h"
#include "common.h"
#include "draw.h"

#include "navier-stokes/perfs.h"

// =======================================================
// Time parameters =======================================

const double tEnd = 1e-2;
const double tStep = 1e-3;

// =======================================================
// Space parameters ======================================

const double jetThickness = 0.16e-3; // 5 or 10mm
const double domainLength = 100. * jetThickness;

// =======================================================
// Physical parameters ===================================
// L for liquid ; G for gaz

scalar f0[], U0[];
char name[80];

// Density
const double rhoL = 1000., rhoG = 1.225;
// Viscosity
const double muL= 1e-3, muG = 1.8e-5;


// Gravity.
const double gravity = 9.81;
// Surface tension
const double sigma = 72e-3;


double Re;
const double We = 1.9; // fixed Weber number to find speed with jetThickness (We = 0.002 -> 6k, smooth around = ???)

double u0;
double lambda; // 2.7mm

// =======================================================
// Mesh parameters =======================================

const int maxlevel = 9;
const double uemax = 1e-3;

// =======================================================
// Boundary conditions ===================================

// Inflow

u.n[left] = dirichlet(f0[] * u0 * (1 - sq(y/jetThickness)));
u.t[left] = dirichlet(0.);
uf.n[left] = dirichlet(f0[] * u0 * (1 - sq(y/jetThickness)));
uf.t[left] = dirichlet(0.);
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

    u0 = sqrt((We * sigma) / (rhoL * jetThickness)); // with We
    //u0 = (Re * muL) / (rhoL * jetThickness); // with Re

    Re = rhoL * u0 * jetThickness / muL;
    lambda = sqrt(sigma / ((rhoL-rhoG)*gravity)) * 1e3; // 2.7mm

    run();
}

// =======================================================
// Initial conditions ====================================

event init(t = 0)
{
    TOLERANCE = 1e-5 [*];

    G.x = gravity;

    fraction(f0, sq(jetThickness) - sq(y));
    f0.refine = f0.prolongation = fraction_refine;
    restriction({f0});

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
        fprintf(stderr, "Re %.1f | We %.2f | u0 %.3f | g %.2f | l %.2f | D %.3f | s %1.1e \n", Re, We, u0, gravity, lambda, jetThickness*1000, sigma);
        fprintf(stderr, "-----------------------------------------------------\n");
        fprintf(stderr, " t / t_max | dt | p.i | f.i | u.i | ncells | perf \n");
        fprintf(stderr, "-----------------------------------------------------\n");

    }

    fprintf(stderr, "%.5f / %.5f | %.2e |  %2d |  %2d |  %2d | %6ld | %.2f\n", t, tEnd, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t);
}

event logEnd(t=tEnd)
{
    fprintf(stderr, "-----------------------------------------------------\n");
    fprintf(stderr, "Re %.1f | We %.2f | u0 %.3f | g %.2f | l %.2f | D %.3f | s %1.1e \n", Re, We, u0, gravity, lambda, jetThickness*1000, sigma);
    fprintf(stderr, "-----------------------------------------------------\n");
}

// =======================================================
// Movie =================================================

event movie(t += tStep; t <= tEnd)
{ 
    view(fov = 22., quat = {0., 0., cos(-pi/4.), cos(pi/4.)}, tx = 0., ty = 0.5, width = 750, height = 750);

    // Figure for f
    clear();
    box();
    draw_vof("f");
    draw_vof ("f", filled = 1, fc = {0.1, 0.1, 0.8});
    sprintf(name, "movie_f_%1.1e.mp4", jetThickness); // TODO : add We
    save(name);

    // Figure for f with cells
    clear();
    cells();
    draw_vof("f");
    draw_vof ("f", filled = 1, fc = {0.1, 0.1, 0.8});
    sprintf(name, "movie_fc_%1.1e.mp4", jetThickness);
    save(name);

    // Figure for u.x
    clear();
    draw_vof("f");
    squares("u.x", linear = true);
    sprintf(name, "movie_u_%1.1e.mp4", jetThickness);
    save(name);

    // Figure for the vorticity
    scalar omega[];
    vorticity(u, omega);
    clear();
    draw_vof("f");
    squares("omega", linear = true);
    sprintf(name, "movie_w_%1.1e.mp4", jetThickness);
    save(name);
}

// =======================================================
// Mesh adaptation =======================================

event adapt(i++)
{
    adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax, uemax}, maxlevel, 3);
}
