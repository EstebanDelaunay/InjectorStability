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
const double tEnd = 0.1e0;
const double tStep = 1e-1;


// =======================================================
// Space parameters ======================================
const double jetThickness = 1e0;
const double domainLength = 50. * jetThickness;


// =======================================================
// Physical parameters ===================================
// L means liquid ; G means gaz

scalar f0[];

// Density
const double rhoL = 1;
const double rhoG = 1e-3;
// Viscosity
const double muL= 1e-3;
const double muG = 1e-5;
// Velocity (defined by Re)
const double u0 = 1e-1;
// Gravity (defined by Fr)
const double gravity = 10;
// Surface tension (defined by We)
const double sigma = 72.5e-5;

const double Re = rhoL * u0 * jetThickness / muL;
const double Fr = u0 * u0 / (gravity * jetThickness);
const double We = rhoL * u0 * u0 * jetThickness / sigma;


// =======================================================
// Mesh parameters =======================================
const int maxlevel = 9; // default maximum level of refinement
const double uemax = 1e-3; // error threshold on velocity


// =======================================================
// Boundary conditions ===================================
// inflow
u.n[left] = dirichlet(f0[] * u0);
u.t[left] = dirichlet(0);
p[left] = neumann(neumann_pressure(0));
f[left] = dirichlet(f0[]);

// Outflow uses standard Neumann/Dirichlet conditions.
u.n[right] = neumann(0);
u.t[right] = neumann(0);
p[right]    = dirichlet(0.);
pf[right]   = dirichlet(0.);
f[right] = neumann(0);


// =======================================================
// main ==================================================
int main()
{
    //periodic(front);

    //init_grid(1 << maxlevel);
    init_grid(64);

    origin(0., -domainLength/2.);
    size(domainLength);

    rho1 = 1.;
    rho2 = rho1 * 1e-3;

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
        fprintf(stderr, "-----------------------------------------------------\n");
        fprintf(stderr, " Re |  Fr |  We |    u0   | gravity |  sigma  \n");
        fprintf(stderr, "%3.f | %3.f | %3.f | %1.1e | %1.1e | %1.1e    \n", Re, Fr, We, u0, gravity, sigma);
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
    view(quat = {0., 0., -cos(pi/4.), cos(pi/4.)}, tx = 0., ty = -0.5);

    clear();
    draw_vof("f");
    squares("f", min=0., max=1.);
    //squares("u.x", min=-2.5, max=2.5);
    //squares("u.x", linear = true);
    //cells();
    //box();
    draw_string(str = "Re", pos = 3, size = 25, lc = {255.,255.,255.}, lw = 4);
    save("movie.mp4");
}


// =======================================================
// Mesh adaptation =======================================
event adapt(i++)
{
    adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax, uemax}, maxlevel, 3);
}
