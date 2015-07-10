#include "grid.hpp"

#include <stdlib.h>

#include "globals.hpp"
#include "messages.hpp"

/* Interaction between tip and surface (interpolation from a grid) */
void interactTipSurfaceFromGrid(void) {

    IVECTOR gp, cgp;
    VECTOR dxyz;
    IVECTOR ai;
    int ix, iy, iz, n;
    double zoffset;
    double *fgr;
    VECTOR C[8];
    double c00, c10, c01, c11, c0, c1, c;

    /* Zero the forces and energy */
    TipSurf_energy = 0.0;
    TipSurf_force = NULL_vector;

    /* What is the nearest grid point to the current position of the tip? */
    retrieveGridPoint(&gp,&dxyz);

    /* Find the surrounding grid points as array indices and retrieve the force values */
    fgr = ForceGridRigid;
    for (ix=0, n=0; ix<=1; ++ix) {
        cgp.x = gp.x + ix;
        for (iy=0; iy<=1; ++iy) {
            cgp.y = gp.y + iy;
            for (iz=0; iz<=1; ++iz) {
                cgp.z = gp.z + iz;
                /* Get the array index */
                retrieveArrayIndex(cgp,&ai);
                /* Retrieve and store the force */
                C[n].x = fgr[ai.x];
                C[n].y = fgr[ai.y];
                C[n].z = fgr[ai.z];
                n++;
            }
        }
    }

    /* This is the order of the coefficients in the C array */
    /*      000 - 001 - 010 - 011 - 100 - 101 - 110 - 111           */

    /* The trilinear interpolation below is based on:                */
    /*  http://en.wikipedia.org/wiki/Trilinear_interpolation */

    /* Construct the force on the tip (in x) */
    c00 = C[0].x*(1-dxyz.x) + C[4].x*(dxyz.x);
    c01 = C[1].x*(1-dxyz.x) + C[5].x*(dxyz.x);
    c10 = C[2].x*(1-dxyz.x) + C[6].x*(dxyz.x);
    c11 = C[3].x*(1-dxyz.x) + C[7].x*(dxyz.x);
    c0 = c00*(1-dxyz.y) + c10*(dxyz.y);
    c1 = c01*(1-dxyz.y) + c11*(dxyz.y);
    TipSurf_force.x = c0*(1-dxyz.z) + c1*(dxyz.z);

    /* Construct the force on the tip (in y) */
    c00 = C[0].y*(1-dxyz.x) + C[4].y*(dxyz.x);
    c01 = C[1].y*(1-dxyz.x) + C[5].y*(dxyz.x);
    c10 = C[2].y*(1-dxyz.x) + C[6].y*(dxyz.x);
    c11 = C[3].y*(1-dxyz.x) + C[7].y*(dxyz.x);
    c0 = c00*(1-dxyz.y) + c10*(dxyz.y);
    c1 = c01*(1-dxyz.y) + c11*(dxyz.y);
    TipSurf_force.y = c0*(1-dxyz.z) + c1*(dxyz.z);

    /* Construct the force on the tip (in z) */
    c00 = C[0].z*(1-dxyz.x) + C[4].z*(dxyz.x);
    c01 = C[1].z*(1-dxyz.x) + C[5].z*(dxyz.x);
    c10 = C[2].z*(1-dxyz.x) + C[6].z*(dxyz.x);
    c11 = C[3].z*(1-dxyz.x) + C[7].z*(dxyz.x);
    c0 = c00*(1-dxyz.y) + c10*(dxyz.y);
    c1 = c01*(1-dxyz.y) + c11*(dxyz.y);
    TipSurf_force.z = c0*(1-dxyz.z) + c1*(dxyz.z);

    /* Go home */
    return;
}

/* Build a 3D force grid (for efficient look ups in case there is no flexibility) */
void build3DForceGrid(void) {

    double mingridspacing, smallestgridspacing;
    double *localforcegrid;
    int ix, iy, iz;
    int i, n, onproc;
    int nx, ny, nz;
    double x, y, z;
    double zoffset;
    FILE *fp;

    /* Check the grid spacing, and adjust if necessary */
    mingridspacing = 9e99;
    smallestgridspacing = 0.01; /* A grid smaller than this (in Angstrom) is silly */
    if (Options.dx < mingridspacing) {
        mingridspacing = Options.dx;
    }
    if (Options.dy < mingridspacing) {
        mingridspacing = Options.dy;
    }
    if (Options.dz < mingridspacing) {
        mingridspacing = Options.dz;
    }
    if (mingridspacing < smallestgridspacing) {
        mingridspacing = smallestgridspacing;
    }
    GridSpacing = mingridspacing;

    /* How many points in the grid in each direction */
    Ngrid.x = (int)((Npoints.x*Options.dx)/GridSpacing);
    Ngrid.y = (int)((Npoints.y*Options.dy)/GridSpacing);
    Ngrid.z = (int)((Npoints.z*Options.dz)/GridSpacing);
    Ngridpoints = (Ngrid.x+1)*(Ngrid.y+1)*(Ngrid.z+1);

    /* Allocate memory for force grids (locally and globally) */
    ForceGridRigid = (double *)malloc(3*Ngridpoints*sizeof(double));
    localforcegrid = (double *)malloc(3*Ngridpoints*sizeof(double));
    for (i=0; i<(3*Ngridpoints); ++i) {
        ForceGridRigid[i] = 0.0;
        localforcegrid[i] = 0.0;
    }

    /* Loop all possible grid points and locally compute the VdW force */
    debugline(RootProc,"Computing 3D force grid (%d grid points)",Ngridpoints);
    zoffset = Options.zhigh - DummyParams.rmin - Ngrid.z*GridSpacing;
    /* First over x */
    for (ix=0; ix<=Ngrid.x; ++ix) {
        /* Where is x */
        Tip_pos.x = ix*GridSpacing;
        /* Then over y */
        for (iy=0; iy<=Ngrid.y; ++iy) {
            /* Where is y */
            Tip_pos.y = iy*GridSpacing;
            /* Which processor are we on? Do we skip this point? */
            n = ix*(Ngrid.y+1) + iy;
            onproc = n % NProcessors;
            if (onproc != Me) {
                continue;
            }
            /* Finally loop over z */
            for (iz=0; iz<=Ngrid.z; ++iz) {
                /* Where is z */
                Tip_pos.z = iz*GridSpacing + zoffset;
                /* Compute the VdW/Coulomb interaction */
                /* We cannot use the function pointer here, this has to be explicit! */
                interactTipSurfaceDirectly();
                /* Which point in the force grid are we storing */
                nx = ((ix*(Ngrid.y+1)) + iy)*(Ngrid.z+1) + iz;
                ny = nx + Ngridpoints;
                nz = ny + Ngridpoints;
                /* Store the force */
                localforcegrid[nx] = TipSurf_force.x;
                localforcegrid[ny] = TipSurf_force.y;
                localforcegrid[nz] = TipSurf_force.z;
            } /* End loop in z */
        } /* End loop in y */
    } /* End loop in x */

    /* First wait and then communicate all the data to all processors */
#if !SERIAL
    MPI_Barrier(Universe);
    MPI_Allreduce(localforcegrid,ForceGridRigid,3*Ngridpoints,MPI_DOUBLE,MPI_SUM,Universe);
#else
    for (i=0; i<(3*Ngridpoints); ++i) {
        ForceGridRigid[i] = localforcegrid[i];
    }
#endif

    /* We probably do not need to write the force grid to file. It will be huge and
         because the grid is computer in parallel, it is relatively fast to recompute it. */

    /* The main processor will write the grid to file */
    //if (Me == RootProc) {
    //  fp = fopen("frcgrid.dat","w");
    //  for (ix=0; ix<=Ngrid.x; ++ix) {
    //      for (iy=0; iy<=Ngrid.y; ++iy) {
    //  for (iz=0; iz<=Ngrid.z; ++iz) {
    //          /* Which points do we need to retrieve? */
    //          nx = ((ix*(Ngrid.y+1)) + iy)*(Ngrid.z+1) + iz;
    //      ny = nx + Ngridpoints;
    //      nz = ny + Ngridpoints;
    //      fprintf(fp,"%d %d %d %10.6f %10.6f %10.6f\n",ix,iy,iz,ForceGridRigid[nx],ForceGridRigid[ny],ForceGridRigid[nz]);
    //  }
    //  }
    //  }
    //  fclose(fp);
    //}

    /* Go home */
    return;
}

/* Retrieve grid point ... */
void retrieveGridPoint(IVECTOR *gp,VECTOR *dxyz) {
    double zoffset;

    /* What is the nearest grid point to the current position of the tip */
    zoffset = Options.zhigh - DummyParams.rmin - Ngrid.z*GridSpacing;
    gp->x = (int)(Tip_pos.x/GridSpacing);
    gp->y = (int)(Tip_pos.y/GridSpacing);
    gp->z = (int)((Tip_pos.z-zoffset)/GridSpacing);

    /* Fix possible out of bounds (only happens at the edges) */
    if (Tip_pos.x<0) {
        gp->x = 0;
    }
    if (Tip_pos.y<0) {
        gp->y = 0;
    }
    if (Tip_pos.z<zoffset) {
        gp->z = 0;
    }
    if (Tip_pos.x>Box.x) {
        gp->x = Ngrid.x-1;
    }
    if (Tip_pos.y>Box.y) {
        gp->y = Ngrid.y-1;
    }
    if (Tip_pos.z>(zoffset+Ngrid.z*GridSpacing)) {
        gp->z = Ngrid.z-1;
    }

    /* And how far are we inside the cube (fractional coordinates) */
    dxyz->x = (Tip_pos.x - ((gp->x)*GridSpacing)) / GridSpacing;
    dxyz->y = (Tip_pos.y - ((gp->y)*GridSpacing)) / GridSpacing;
    dxyz->z = ((Tip_pos.z-zoffset) - ((gp->z)*GridSpacing)) / GridSpacing;

    /* Fix possible out of bounds (only happens at the edges) */
    if (Tip_pos.x<0) {
        dxyz->x = 0.0;
    }
    if (Tip_pos.y<0) {
        dxyz->y = 0.0;
    }
    if (Tip_pos.z<zoffset) {
        dxyz->z = 0.0;
    }
    if (Tip_pos.x>Box.x) {
        dxyz->x = 1.0;
    }
    if (Tip_pos.y>Box.y) {
        dxyz->y = 1.0;
    }
    if (Tip_pos.z>(zoffset+Ngrid.z*GridSpacing)) {
        dxyz->z = 1.0;
    }

    /* Go home */
    return;
}

/* ... and corresponding array indices */
void retrieveArrayIndex(IVECTOR gp, IVECTOR *ai) {

    /* Compute the array index based on the given grid point */
    ai->x = ((gp.x*(Ngrid.y+1)) + gp.y)*(Ngrid.z+1) + gp.z;
    ai->y = ai->x + Ngridpoints;
    ai->z = ai->y + Ngridpoints;

    /* Go home */
    return;
}
