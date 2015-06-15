#include "physics.hpp"

#include <math.h>

#include "globals.hpp"
#include "grid.hpp"


void moveTip(void) {

    int n, nmax, check;
    int i, ix, iy, iz;
    double x, y, z;
    double angle, minangle, maxforce, dd;
    double e, ediff, eold, fnorm;
    VECTOR f, d;
    VECTOR *ftip;
    int nxy, onproc, bufsize, nsr;
    BUFFER *sendbuf, *recvbuf;
    double checkperc, curperc;

    FILE *tfp;
    char dump[NAME_LENGTH];

    // [> Create a force vector (for real time analysis) <]
    ftip = (VECTOR *)malloc((Npoints.z+1)*sizeof(VECTOR));

    // [> DEBUG DEBUG DEBUG DEBUG DEBUG <]
    //Npoints.x = Npoints.y = 10;

    // [> Some initialization <]
    Ntotal = 0;
    checkperc = curperc = 0.10;
    debugline(simulation, RootProc,"Simulation run started now");

    // [> Initialize the storage send and receive buffers <]
    nxy = (Npoints.x + 1)*(Npoints.y + 1);
    bufsize = Options.bufsize * (Npoints.z+1);
    sendbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    recvbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    nsr = 0;

    // [> Loop x <]
    for (ix=0; ix<=Npoints.x; ++ix) {
        x = ix*Options.dx; [> Current x <]

        // [> Loop y <]
        for (iy=0; iy<=Npoints.y; ++iy) {
            y = iy*Options.dy; [> Current y <]

            // [> Check the progress and report every so often <]
            n = ix*(Npoints.y+1) + iy;
            if ( (Me == RootProc) && ((((double)n)/(nxy)) >= curperc) ) {
                debugline(simulation, RootProc,"Finished approximately %4.1f %% of the simulation",100*curperc);
                curperc += checkperc;
            }

            // [> Compute on which processor this x,y combination should be run <]
            onproc = n % NProcessors;
            if (onproc != Me) {
                continue;
            }
            PointsOnProc[Me]++;

            // [> Position the tip far above the surface <]
            Dummy_pos.x = x;
            Dummy_pos.y = y;
            Dummy_pos.z = Options.zhigh + Options.dz; [> Plus dz to allow for initial subtraction <]
            Tip_pos.x = Dummy_pos.x;
            Tip_pos.y = Dummy_pos.y;
            Tip_pos.z = Dummy_pos.z - DummyParams.rmin;

            // [> If the molecule is flexible, reset it to its original position, before beginning the approach <]
            if (Options.flexible) {
                for (i=0; i<Natoms; ++i) {
                    Surf_pos[i].x = Surf_pos_org[i].x;
                    Surf_pos[i].y = Surf_pos_org[i].y;
                    Surf_pos[i].z = Surf_pos_org[i].z;
                }
            }

            // [> Approach and optimize <]
            nmax = 0;
            minangle = 9e99;
            maxforce = -9e99;
            for (iz=0; iz<=Npoints.z; ++iz) {
                z = Options.zhigh - iz*Options.dz; [> Current z <]

                // [> Move tip and dummy atom toward the surface <]
                Dummy_pos.z -= Options.dz;
                Tip_pos.z -= Options.dz;

                // [> Collect the force <]
                ftip[iz] = NULL_vector;

                // [> Relax/Minimize the configuration <]
                ediff = 5*Options.etol;
                for (n=0; n<Options.maxsteps; ++n) {
                    nmax++;

                    // [> Compute all interaction energies <]
                    interactTipSurface();
                    interactTipDummy();
                    interactTipHarmonic();

                    // [> Total energy and force computed <]
                    e = TipSurf_energy + TipDummy_energy + TipHarmonic_energy;
                    f.x = TipSurf_force.x + TipDummy_force.x + TipHarmonic_force.x;
                    f.y = TipSurf_force.y + TipDummy_force.y + TipHarmonic_force.y;
                    f.z = TipSurf_force.z + TipDummy_force.z + TipHarmonic_force.z;
                    fnorm = sqrt(f.x*f.x + f.y*f.y + f.z*f.z);

                    // [> Energy difference <]
                    if (n>0) {
                        ediff = e - eold;
                    }
                    eold = e;

                    // [> Are the forces/energies tolerable <]
                    if (Options.minterm == MIN_E) {
                        check = (fabs(ediff)<Options.etol);
                    }
                    else if (Options.minterm == MIN_F) {
                        check = (fabs(fnorm)<Options.ftol);
                    }
                    else if (Options.minterm == MIN_EF) {
                        check = ((fabs(ediff)<Options.etol)&&(fabs(fnorm)<Options.ftol));
                    }
                    if (check) {
                        break;
                    }

                    // [> If they are not tolerable, update position of the tip atom based on the force <]
                    Tip_pos.x += Options.cfac * f.x;
                    Tip_pos.y += Options.cfac * f.y;
                    Tip_pos.z += Options.cfac * f.z;

                    // [> If the molecule is supposed to behave flexible, we need to update its positions too <]
                    if (Options.flexible) {
                        updateFlexibleMolecule();
                    }

                } // [> End minimization loop <]

                // [> Compute some other interesting data <]
                d.x = Tip_pos.x - x;
                d.y = Tip_pos.y - y;
                d.z = Tip_pos.z - z;
                dd = sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
                angle = atan2(sqrt(d.x*d.x + d.y*d.y),d.z)*(180.0/PI);
                if (angle<minangle) {
                    minangle = angle;
                }
                if (fnorm>maxforce) {
                    maxforce = fnorm;
                }

                // [> Collect the final forces on the tip because of the surface <]
                ftip[iz].x = TipSurf_force.x;
                ftip[iz].y = TipSurf_force.y;
                ftip[iz].z = TipSurf_force.z;

                // [> Store data in send buffers <]
                sendbuf[nsr].ix = ix;
                sendbuf[nsr].iy = iy;
                sendbuf[nsr].iz = iz;
                sendbuf[nsr].n  = n;
                sendbuf[nsr].pos.x = x;
                sendbuf[nsr].pos.y = y;
                sendbuf[nsr].pos.z = z;
                sendbuf[nsr].f.x = TipSurf_force.x;
                sendbuf[nsr].f.y = TipSurf_force.y;
                sendbuf[nsr].f.z = TipSurf_force.z;
                sendbuf[nsr].d.x = d.x;
                sendbuf[nsr].d.y = d.y;
                sendbuf[nsr].d.z = d.z;
                sendbuf[nsr].dd = dd;
                sendbuf[nsr].e = TipSurf_energy;
                sendbuf[nsr].angle = angle;
                nsr++;

            } // [> End loop in z <]

            // [> Dump to file (only when the buffer is full, this can and should only happen after a full z approach) <]
            if (nsr==bufsize) {
                dumpToFiles(sendbuf,recvbuf,bufsize);
                nsr = 0;
            }

            // [> Compute the frequency shift for the given F(z) <]
            //computeDeltaF(x,y,ftip);

            if ((ix==(Npoints.x/2)) && (iy==(Npoints.y/2))) {
                sprintf(dump,"test-%d-%d.xyz",ix,iy);
                tfp = fopen(dump,"w");
                fprintf(tfp,"%d\n\n",Natoms);
                for (i=0; i<Natoms; ++i) {
                    fprintf(tfp,"%s %8.4f %8.4f %8.4f\n",Surf_type[i],Surf_pos[i].x,Surf_pos[i].y,Surf_pos[i].z);
                }
                fclose(tfp);
            }

            // [> Keep track of counting <]
            Ntotal += nmax;

        } // [> End loop in y <]

    } // [> End loop in x <]

    // [> Dump to file (if it happened that the buffer contains anything at all <]
#if !SERIAL
    MPI_Allreduce(&nsr,&n,1,MPI_INT,MPI_SUM,Universe);
#else
    n = nsr;
#endif
    if (n>0) {
        dumpToFiles(sendbuf,recvbuf,nsr);
    }

    // [> Say one last thing <]
    debugline(simulation, RootProc,"Finished %5.1f %% of the simulation",100*curperc);

    // [> Go home <]
    return;
}

/* Interaction between tip and surface */
void interactTipSurfaceDirectly(void) {

    int i;
    double dx, dy, dz, rsqt, dr, dexp, frcfac;
    double e, fx, fy, fz;
    double epair, fpair, terma, termb, termc, sr6, sr12;

    /* Zero the forces and energy */
    TipSurf_energy = 0.0;
    TipSurf_force = NULL_vector;

    /* Zero the forces for the molecule if it's supposed to be flexible */
    if (Options.flexible) {
        for (i=0; i<Natoms; ++i) {
            Surf_force[i].x = 0.0;
            Surf_force[i].y = 0.0;
            Surf_force[i].z = 0.0;
        }
    }

    /* Loop all surface particles */
    for (i=0; i<Natoms; ++i) {
        /* Compute distance (components) */
        dx = Tip_pos.x - Surf_pos[i].x;
        dy = Tip_pos.y - Surf_pos[i].y;
        dz = Tip_pos.z - Surf_pos[i].z;
        rsqt = dx*dx + dy*dy + dz*dz;
        /* Compute the Van der Waals contribution factors */
        if (TipSurfParams[i].morse) {
            /* The LJ is replaced by a Morse potential */
            dr = sqrt(rsqt) - TipSurfParams[i].re;
            dexp = exp(-TipSurfParams[i].a * dr);
            frcfac = 2.0 * TipSurfParams[i].De * TipSurfParams[i].a;
            epair = TipSurfParams[i].De * ( dexp*dexp - 2.0*dexp + 1 );
            fpair = frcfac * ( dexp*dexp - dexp )/sqrt(rsqt);
        }
        else {
            /* The Lennard-Jones 12-6 interaction coefficients */
            sr6 = rsqt*rsqt*rsqt;
            sr12 = sr6*sr6;
            terma = TipSurfParams[i].es12/sr12;
            termb = TipSurfParams[i].es6/sr6;
            epair = terma-termb;
            fpair = (12*terma-6*termb)/rsqt;
        }
        /* The Coulomb interaction coefficients (only if coulomb is on) */
        if (Options.coulomb) {
            termc = TipSurfParams[i].qq/sqrt(rsqt);
        }
        else {
            termc = 0.0;
        }
        /* The interaction energy */
        TipSurf_energy += epair + termc;
        /* The interaction force */
        TipSurf_force.x += (fpair+(termc/rsqt))*dx;
        TipSurf_force.y += (fpair+(termc/rsqt))*dy;
        TipSurf_force.z += (fpair+(termc/rsqt))*dz;
        /* If the molecule is flexible, also store its forces */
        if (Options.flexible) {
            Surf_force[i].x -= (fpair+(termc/rsqt))*dx;
            Surf_force[i].y -= (fpair+(termc/rsqt))*dy;
            Surf_force[i].z -= (fpair+(termc/rsqt))*dz;
        }
    }

    /* Go home */
    return;
}

/* Interaction between tip and dummy atom */
void interactTipDummy(void) {

    double dx, dy, dz, rsqt, dr, dexp, frcfac;
    double epair, fpair, terma, termb, termc, sr6, sr12;

    /* Compute distance (components) */
    dx = Tip_pos.x - Dummy_pos.x;
    dy = Tip_pos.y - Dummy_pos.y;
    dz = Tip_pos.z - Dummy_pos.z;
    rsqt = dx*dx + dy*dy + dz*dz;
    /* The Lennard-Jones interaction coefficients */
    sr6 = rsqt*rsqt*rsqt;
    sr12 = sr6*sr6;
    terma = DummyParams.es12/sr12;
    termb = DummyParams.es6/sr6;
    epair = terma-termb;
    fpair = (12*terma-6*termb)/rsqt;
    /* Now check if the LJ is replaced by a Morse potential */
    if (DummyParams.morse) {
        dr = sqrt(rsqt) - DummyParams.re;
        dexp = exp(-DummyParams.a * dr);
        frcfac = 2.0 * DummyParams.De * DummyParams.a;
        epair = DummyParams.De * ( dexp*dexp - 2.0*dexp + 1 );
        fpair = frcfac * ( dexp*dexp - dexp )/sqrt(rsqt);
    }
    /* The interaction energy */
    TipDummy_energy = epair;
    /* The interaction force */
    TipDummy_force.x = fpair*dx;
    TipDummy_force.y = fpair*dy;
    TipDummy_force.z = fpair*dz;

    /* Go home */
    return;
}

/* Interaction between tip and harmonic constraint */
void interactTipHarmonic(void) {

    double dx, dy, r, dr, rk, fharm;

    /* Compute distance (components, but not in z!) */
    dx = Tip_pos.x - Dummy_pos.x;
    dy = Tip_pos.y - Dummy_pos.y;
    r = sqrt(dx*dx + dy*dy);
    /* Compute the harmonic coefficients */
    dr = r - Harmonic.r0;
    rk = Harmonic.k*dr;
    /* The interaction energy */
    TipHarmonic_energy = rk*dr;
    /* The interaction force */
    if (r>TOLERANCE) {
        fharm = (-2*rk/r);
    }
    else {
        fharm = 0.0;
    }
    TipHarmonic_force.x = fharm*dx;
    TipHarmonic_force.y = fharm*dy;
    TipHarmonic_force.z = 0.0;

    /* Go home */
    return;
}

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

/* void computeDeltaF(double x, double y, VECTOR *ftip) { */

/*   int i; */
/*   double z; */

/*   /\* Wait! *\/ */
/*   MPI_Barrier(Universe); */

/*   /\* Copy data *\/ */
/*   for (i=0; i<=Npoints.z; ++i) { */
/*       z = Options.zhigh - i*Options.dz; */
/*       fprintf(stdout,"@ %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",x,y,z,ftip[i].x,ftip[i].y,ftip[i].z); */
/*   } */

/*   /\* Go home *\/ */
/*   return; */
/* } */
