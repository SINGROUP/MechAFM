/****************************
 **                        **
 **  Mechanical AFM Model  **
 **                        **
 **    Based on paper:     **
 **   Hapala et al, PRB    **
 **    90:085421 (2014)    **
 **                        **
 **  C+MPI Implementation  **
 **       2014/2015        **
 **      Peter Spijker     **
 **    Aalto University    **
 ** peter.spijker@aalto.fi **
 **                        **
 ****************************/

/****
TO DO:
- instead of fixing atoms to mimic support, let's put a wall potential which prevents atoms moving too far "down"
- VDW-E interactions between molecules (not intramolecular), needs molecule identification and clustering
- implement better minimization algorithm
- create DF image on the fly
- implement the Hartree potential (grid based)
****/

/* Load system headers */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <glob.h>
#include <sys/time.h>
#if !SERIAL
    #include <mpi.h>
#endif

#include "flexible.hpp"
#include "globals.hpp"
#include "grid.hpp"
#include "messages.hpp"
#include "utility.hpp"
#include "parse.hpp"
#include "physics.hpp"
#include "simulation.hpp"
#include "system.hpp"

/**********************
 ** GLOBAL VARIABLES **
 **********************/

/* Function pointers */
void (*interactTipSurface)(void);       // For the tip surface interaction
void interactTipSurfaceDirectly(void);
void interactTipSurfaceFromGrid(void);

/*************************************************
 ** EVERYTHING TO DO WITH THE SIMULATION ITSELF **
 *************************************************/

/* Set some global variables we need and open the file streams */
void openUniverse(void) {

    int i, n;
    double z;
    char outfile[NAME_LENGTH];

    /* How many points to compute the tip relaxation on */
    Npoints.x = (int)(Box.x/Options.dx);
    Npoints.y = (int)(Box.y/Options.dy);
    Npoints.z = (int)((Options.zhigh-Options.zlow)/Options.dz);
    n = (Npoints.x + 1) * (Npoints.y + 1) * (Npoints.z + 1);
    debugline(RootProc,"3D data grid is: %d x %d x %d (%d in total)",1+Npoints.x,1+Npoints.y,1+Npoints.z,n);

#if !SERIAL
    /* Wait! */
    MPI_Barrier(Universe);
#endif

    /* If we which to precompute the force grid, do it now */
    if (Options.rigidgrid) {
        build3DForceGrid();
    }

    /* Open all the file streams (one for every z point) [ONLY ON ROOT PROCESSOR] */
    if (Me == RootProc) {
        FStreams = (FILE **)malloc((Npoints.z+1)*sizeof(FILE *));
        for (i=0; i<=Npoints.z; ++i) {
            z = Options.zhigh - i*Options.dz;
            if ( Options.gzip == TRUE ) {
                sprintf(outfile,"gzip -6 > scan-%06.3f.dat.gz",z);
                FStreams[i] = popen(outfile,"w");
            }
            else {
                sprintf(outfile,"scan-%06.3f.dat",z);
                FStreams[i] = fopen(outfile,"w");
            }
        }
    }

    /* Note the time */
    gettimeofday(&TimeStart, NULL);

    /* Go home */
    return;
}

/* Now move the tip! */
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

    /* Create a force vector (for real time analysis) */
    ftip = (VECTOR *)malloc((Npoints.z+1)*sizeof(VECTOR));

    /* DEBUG DEBUG DEBUG DEBUG DEBUG */
    //Npoints.x = Npoints.y = 10;

    /* Some initialization */
    Ntotal = 0;
    checkperc = curperc = 0.10;
    debugline(RootProc,"Simulation run started now");

    /* Initialize the storage send and receive buffers */
    nxy = (Npoints.x + 1)*(Npoints.y + 1);
    bufsize = Options.bufsize * (Npoints.z+1);
    sendbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    recvbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    nsr = 0;

    /* Loop x */
    for (ix=0; ix<=Npoints.x; ++ix) {
        x = ix*Options.dx; /* Current x */

        /* Loop y */
        for (iy=0; iy<=Npoints.y; ++iy) {
            y = iy*Options.dy; /* Current y */

            /* Check the progress and report every so often */
            n = ix*(Npoints.y+1) + iy;
            if ( (Me == RootProc) && ((((double)n)/(nxy)) >= curperc) ) {
                debugline(RootProc,"Finished approximately %4.1f %% of the simulation",100*curperc);
                curperc += checkperc;
            }

            /* Compute on which processor this x,y combination should be run */
            onproc = n % NProcessors;
            if (onproc != Me) {
                continue;
            }
            PointsOnProc[Me]++;

            /* Position the tip far above the surface */
            Dummy_pos.x = x;
            Dummy_pos.y = y;
            Dummy_pos.z = Options.zhigh + Options.dz; /* Plus dz to allow for initial subtraction */
            Tip_pos.x = Dummy_pos.x;
            Tip_pos.y = Dummy_pos.y;
            Tip_pos.z = Dummy_pos.z - DummyParams.rmin;

            /* If the molecule is flexible, reset it to its original position, before beginning the approach */
            if (Options.flexible) {
                for (i=0; i<Natoms; ++i) {
                    Surf_pos[i].x = Surf_pos_org[i].x;
                    Surf_pos[i].y = Surf_pos_org[i].y;
                    Surf_pos[i].z = Surf_pos_org[i].z;
                }
            }

            /* Approach and optimize */
            nmax = 0;
            minangle = 9e99;
            maxforce = -9e99;
            for (iz=0; iz<=Npoints.z; ++iz) {
                z = Options.zhigh - iz*Options.dz; /* Current z */

                /* Move tip and dummy atom toward the surface */
                Dummy_pos.z -= Options.dz;
                Tip_pos.z -= Options.dz;

                /* Collect the force */
                ftip[iz] = NULL_vector;

                /* Relax/Minimize the configuration */
                ediff = 5*Options.etol;
                for (n=0; n<Options.maxsteps; ++n) {
                    nmax++;

                    /* Compute all interaction energies */
                    interactTipSurface();
                    interactTipDummy();
                    interactTipHarmonic();

                    /* Total energy and force computed */
                    e = TipSurf_energy + TipDummy_energy + TipHarmonic_energy;
                    f.x = TipSurf_force.x + TipDummy_force.x + TipHarmonic_force.x;
                    f.y = TipSurf_force.y + TipDummy_force.y + TipHarmonic_force.y;
                    f.z = TipSurf_force.z + TipDummy_force.z + TipHarmonic_force.z;
                    fnorm = sqrt(f.x*f.x + f.y*f.y + f.z*f.z);

                    /* Energy difference */
                    if (n>0) {
                        ediff = e - eold;
                    }
                    eold = e;

                    /* Are the forces/energies tolerable */
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

                    /* If they are not tolerable, update position of the tip atom based on the force */
                    Tip_pos.x += Options.cfac * f.x;
                    Tip_pos.y += Options.cfac * f.y;
                    Tip_pos.z += Options.cfac * f.z;

                    /* If the molecule is supposed to behave flexible, we need to update its positions too */
                    if (Options.flexible) {
                        updateFlexibleMolecule();
                    }

                } /* End minimization loop */

                /* Compute some other interesting data */
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

                /* Collect the final forces on the tip because of the surface */
                ftip[iz].x = TipSurf_force.x;
                ftip[iz].y = TipSurf_force.y;
                ftip[iz].z = TipSurf_force.z;

                /* Store data in send buffers */
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

            } /* End loop in z */

            /* Dump to file (only when the buffer is full, this can and should only happen after a full z approach) */
            if (nsr==bufsize) {
                dumpToFiles(sendbuf,recvbuf,bufsize);
                nsr = 0;
            }

            /* Compute the frequency shift for the given F(z) */
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

            /* Keep track of counting */
            Ntotal += nmax;

        } /* End loop in y */

    } /* End loop in x */

    /* Dump to file (if it happened that the buffer contains anything at all */
#if !SERIAL
    MPI_Allreduce(&nsr,&n,1,MPI_INT,MPI_SUM,Universe);
#else
    n = nsr;
#endif
    if (n>0) {
        dumpToFiles(sendbuf,recvbuf,nsr);
    }

    /* Say one last thing */
    debugline(RootProc,"Finished %5.1f %% of the simulation",100*curperc);

    /* Go home */
    return;
}

/* Close all the file streams */
void closeUniverse(void) {

    int i;

#if !SERIAL
    /* Wait! */
    MPI_Barrier(Universe);
#endif

    /* Close each separate file stream [ONLY ON ROOT PROCESSOR] */
    if (Me == RootProc) {
        for (i=0; i<=Npoints.z; ++i) {
            if ( Options.gzip == TRUE ) {
                pclose(FStreams[i]);
            }
            else {
                fclose(FStreams[i]);
            }
        }
    }

    /* Go home */
    return;
}

/********************
 ** FINAL THOUGHTS **
 ********************/

void finalize(void) {

    int n, nsum;
    double dtime, timesum;

    /* Note the time */
    gettimeofday(&TimeEnd, NULL);

    /* Time difference */
    dtime = (TimeEnd.tv_sec - TimeStart.tv_sec + (TimeEnd.tv_usec - TimeStart.tv_usec)/1e6);
    timesum = 0.0;
#if !SERIAL
    MPI_Allreduce(&dtime,&timesum,1,MPI_DOUBLE,MPI_SUM,Universe);
#else
    timesum += dtime;
#endif

    /* Collect number of steps from all processors */
    nsum = 0;
#if !SERIAL
    MPI_Allreduce(&Ntotal,&nsum,1,MPI_INT,MPI_SUM,Universe);
#else
    nsum += Ntotal;
#endif

    /* Print some miscelleneous information */
    debugline(RootProc,"Simulation run finished");
    debugline(RootProc,"Statistics:");
    n = (Npoints.x + 1) * (Npoints.y + 1) * (Npoints.z + 1);
    debugline(RootProc,"    Computed %ld tip positions",n);
    debugline(RootProc,"    Needed %ld minimization steps in total",nsum);
    debugline(RootProc,"    Which means approximately %.2f minimization steps per tip position",((double)nsum/n));
    debugline(RootProc,"    The simulation wall time is %.2f seconds",timesum);
    debugline(RootProc,"    The entire simulation took %.2f seconds",dtime);
    debugline(RootProc,"");

    /* Go home */
    return;
}

/***************************
 ** THE PARALLEL UNIVERSE **
 ***************************/

/* Initialize our parallel world */
void openParallelUniverse(int argc, char *argv[], Simulation& simulation) {
#if !SERIAL
    /* Start MPI */
    MPI_Init(&argc,&argv);
#endif
    /* Determine the size of the universe and which processor we are on */
    simulation.root_processor_ = 0;
#if !SERIAL
    simulation.universe = MPI_COMM_WORLD;
    MPI_Comm_rank(Universe, &simulation.me_);
    MPI_Comm_size(Universe, &simulation.n_processors_);
#else
    simulation.me_ = 0;
    simulation.n_processors_ = 1;
#endif
    return;
}

/* Terminate our parallel worlds */
void closeParallelUniverse(void) {

    int i;
    int *pop;

    /* How many x,y points on each processor */
#if !SERIAL
    MPI_Barrier(Universe);
#endif
    pop = (int *)malloc(NProcessors*sizeof(int));
    for (i=0; i<NProcessors; ++i) {
        pop[i] = 0;
    }
#if !SERIAL
    MPI_Allreduce(PointsOnProc,pop,NProcessors,MPI_INT,MPI_SUM,Universe);
#else
    pop[Me] += PointsOnProc[Me];
#endif
    debugline(RootProc,"How many x,y points did each processor handle:");
    for (i=0; i<NProcessors; ++i) {
        debugline(RootProc,"    Processor %2d: %6d x,y points",i,pop[i]);
    }

#if !SERIAL
    /* Close MPI */
    MPI_Finalize();
#endif

    /* Go home */
    return;
}

/**********************
 ** THE MAIN ROUTINE **
 **********************/

int main(int argc, char *argv[]) {

    Simulation simulation;
    /* Set up the parallel routines */
    openParallelUniverse(argc, argv, simulation);

    /* Initialize the simulation */
    parseCommandLine(argc, argv, simulation);        /* Read the command line */
    readInputFile();                            /* Read input file */
    readXYZFile();                              /* Read the XYZ file */
    readParameterFile();                    /* Read the parameter file */

    /* If the molecule is flexible, build the topology */
    if (Options.flexible) {
        buildTopology();
    }

    /* The simulation itself */
    openUniverse();
    moveTip();
    closeUniverse();


    /* Some final thoughts */
    finalize();

    /* And stop the parallel routines properly */
    closeParallelUniverse();

    /* Done */
    return 0;
}
