/**
# Atomisation of a pulsed liquid jet

A dense cylindrical liquid jet is injected into a stagnant lighter
phase (density ratio 1/27.84). The inflow velocity is modulated
sinusoidally to promote the growth of primary shear
instabilities. Surface tension is included and ultimately controls the
characteristic scale of the smallest droplets.

We solve the two-phase Navier--Stokes equations with surface
tension. We need the *tag()* function to count the number of
droplets. We generate animations online using Basilisk View. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"
#include "common.h"

// =======================================================
// Time parameters =======================================
const double tEnd = 5e0;
const double tStep = 5e-1;

// =======================================================
// Space parameters =======================================
const double domainLength = 6e0 [1];


/**
We define the radius of the jet, the initial jet length, the Reynolds
number and the surface tension coefficient. */

double radius = 1. / 12.;
double length = 0.025;
double Re = 5800;
double SIGMA = 3e-5;

/**
The default maximum level of refinement is 10 and the error threshold
on velocity is 0.1. */

int maxlevel = 7;
double uemax = 0.1;

/**
To impose boundary conditions on a disk we use an auxilliary volume 
fraction field *f0* which is one inside the cylinder and zero outside. 

We then set an oscillating inflow velocity on the
left-hand-side and free outflow on the right-hand-side. */

double u0 = 1., au = 0.05, T0 = 0.1;

scalar f0[];
u.n[left] = dirichlet(f0[] * (u0 + au * sin(2. * pi * t / T0)));
u.t[left] = dirichlet(0);


#if dimension > 2
u.r[left] = dirichlet(0);
#endif

p[left] = neumann(0);
f[left] = f0[];

u.n[right] = neumann(0);
p[right] = dirichlet(0);

// =======================================================
// main ==================================================

/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main(int argc, char *argv[])
{
    periodic(top);

    if (argc > 1)
        maxlevel = atoi(argv[1]);
    if (argc > 2)
        uemax = atof(argv[2]);

    /**
    The initial domain is discretised with $64^3$ grid points. We set
    the origin and domain size. */

    init_grid(64);
    origin(0, -domainLength/2., -domainLength/2.);
    size(domainLength);

    /**
    We set the density and viscosity of each phase as well as the
    surface tension coefficient and start the simulation. */

    rho1 = 1. [0];
    rho2 = rho1 / 27.84;

    mu1 = 2. * u0 * radius / Re * rho1;
    mu2 = 2. * u0 * radius / Re * rho2;

    f.sigma = SIGMA;

    run();
}

// =======================================================
// Initial conditions ====================================

event init(t = 0)
{
    if (!restore(file = "restart"))
    {

        /* We use a static refinement down to *maxlevel* in a cylinder 1.2
        times longer than the initial jet and twice the radius. */

        refine(x < 1.2 * length && sq(y) + sq(z) < 2. * sq(radius) && level < maxlevel);

        /* We initialise the auxilliary volume fraction field for a cylinder of
        constant radius. */

        fraction(f0, 1. - );
        f0.refine = f0.prolongation = fraction_refine;
        restriction({f0}); // for boundary conditions on levels

        /* We then use this to define the initial jet and its velocity. */
        foreach ()
        {
            f[] = f0[] * (x < length);
            u.x[] = u0 * f[];
        }
    }
}

// =======================================================
// Log ===================================================

event logfile(i++)
{
    if (i == 0)
    {
        fprintf(stderr, "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
    }

    fprintf(stderr, "%g %g %d %d %d %ld %g %g\n",
            t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t, perf.speed);
}

// =======================================================
// Movie =================================================

event movie(t += tStep; t <= tEnd)
{
    #if dimension == 2
        scalar omega[];
        vorticity(u, omega);
        view(tx = -0.5);
        clear();
        draw_vof("f");
        squares("omega", linear = true, spread = 10);
        box();
    #else  // 3D
        scalar pid[];
        foreach ()
            pid[] = fmod(pid() * (npe() + 37), npe());
        view(camera = "iso",
            fov = 14.5, tx = -0.418, ty = 0.288,
            width = 1600, height = 1200);
        clear();
        draw_vof("f");
    #endif // 3D
    save("movie.mp4");
}

// =======================================================
// Mesh adaptation =======================================

event adapt(i++)
{
    adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax, uemax}, maxlevel);
}