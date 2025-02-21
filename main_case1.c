/* =======================================================================================================

CC99='mpicc -std=c99' qcc -Wall -O2 -autolink -D_MPI=1 mainMD.c -o mainMD -lm -lgsl -lgslcblas -lfb_tiny ; mpirun -np 8 ./mainMD

qcc -autolink -Wall -O2 mainMD.c -o mainMD -lm -lgsl -lcblas -lfb_tiny ; ./mainMD

======================================================================================================= */

#include "grid/octree.h"

#include "navier-stokes/centered.h"

#include "two-phase.h"
#include "tension.h"
#include "reduced.h"

#include "navier-stokes/conserving.h"

#include "view.h"
#include "common.h"
#include "draw.h"

#include "navier-stokes/perfs.h"

#include "signature.h"

#define SIGN_LEV 7

// =======================================================
// Time parameters =======================================

const double tEnd = 0.15e0;
const double tStep = 5e-3;

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
// Gravity.
const double gravity = 9.81;
// Surface tension
const double sigma = 72e-3;
// Weber number
const double We = 10;

double Re, u0, lambda;


// =======================================================
// Mesh parameters =======================================

const int maxlevel = 8;
const double uemax = 1e-3;

scalar phii[], M[];


// =======================================================
// Boundary conditions ===================================

// =======================================================
// Inflow

u.n[left] = dirichlet(f0[] * u0 * (1 - sq(y/jetThickness)));
u.t[left] = dirichlet(0.);
u.r[left] = dirichlet(0.);

uf.n[left] = dirichlet(f0[] * u0 * (1 - sq(y/jetThickness)));
uf.t[left] = dirichlet(0.);
uf.r[left] = dirichlet(0.);

f[left]   = f0[];

// =======================================================
// Outflow

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
u.r[right] = neumann(0.);

p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

f[right]   = dirichlet(0.);

// =======================================================
// main ==================================================

int main()
{
    periodic(front);

    init_grid(1 << 6);

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
    TOLERANCE = 1e-7 [*];

    G.x = gravity;

    fraction(f0, sq(jetThickness) - sq(y));
    f0.refine = f0.prolongation = fraction_refine;
    restriction({f0});

    foreach ()
    {
        f[] = f0[]  * (x < 0.015);
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

event movie (t += tStep)
{
    view(fov = 38, camera = "iso", ty = .4, width = 750, height = 750, bg = {1,1,1}, samples = 4);

    clear();

    squares ("u.y", linear = true);
    squares ("u.x", linear = true, n = {1,0,0});

    // Vorticity
    scalar omega[];
    vorticity (u, omega);
    squares("omega", linear = true, n = {0,1,0});

    // Volume fraction
    draw_vof("f");

    char name[80];
    sprintf(name, "movie_case1_%.3f.png", t);
    save(name);
}

// =======================================================
// Manifold Death Algorithm ==============================


static inline void myrestrict(Point point, scalar s)
{
  s[] = 0;
}

static inline void myprolongation(Point point, scalar s)
{
  foreach_child()
    s[] = 0;
}


scalar thin_structure[];
scalar sign[];

event compute_thin_structure(t += 0.01)
{
    thin_structure.prolongation = myprolongation;
    thin_structure.restriction = myrestrict;
    thin_structure.refine = refine_injection;


    foreach()
    {
        phii[] = 2*f[] - 1;
        sign[] = f[] > 0.5 ? 2 : -1;
    }

    // We choose at which level (`l_sign`) we want to compute the quadratic moments and the signature.
    // Then, the indicator function `phii` is restricted at level `l_sign`.

    int l_sign = SIGN_LEV;

    for (int ilev = depth() - 1; ilev >= l_sign; ilev--)
        foreach_level(ilev){
        if(is_refined(cell))
        restriction_average(point, phii);
        }

    compute_signature_neigh_level (f, phii, sign, l_sign);

    /**
    The signature `sign` is available only at the level `l_sign`.
    We need to prolong it onto the finest grid. */

    for (int ilev = l_sign; ilev < depth(); ilev++)
        foreach_level(ilev)
        {
            sign.prolongation = phii.prolongation = refine_injection;
            if(is_refined(cell))
            {
                sign.prolongation (point, sign);
                phii.prolongation (point, phii);
            }
        }

        foreach()
        if ((f[] > 0.01 && sign[] == -1) || (f[] < 0.99 && sign[] == 2))
            sign[] = 1;

    foreach()
    {
        if(sign[] == 1)
        {
            thin_structure[] = 1;
        }
        else
        {
            thin_structure[] = 0;
        }
    }

    boundary({thin_structure});
    restriction({thin_structure});
}

// =======================================================
// Mesh adaptation =======================================

event adapt(i++)
{
    adapt_wavelet({f, u, thin_structure}, (double[]){uemax, uemax, uemax, uemax, uemax}, maxlevel, 3);
}