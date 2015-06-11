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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#if !SERIAL
    #include <mpi.h>
#endif
#include <chrono>

#include "globals.hpp"
#include "messages.hpp"
#include "parse.hpp"
#include "simulation.hpp"
#include "system.hpp"

using namespace std;

/*************************************************
 ** EVERYTHING TO DO WITH THE SIMULATION ITSELF **
 *************************************************/

/* Set some global variables we need and open the file streams */
void openUniverse(Simulation& simulation) {

    InputOptions& options = simulation.options_;

    int n;
    double z;
    char outfile[NAME_LENGTH];

    /* How many points to compute the tip relaxation on */
    simulation.n_points_.x = (int)(simulation.box_.x / options.dx);
    simulation.n_points_.y = (int)(simulation.box_.y / options.dy);
    simulation.n_points_.z = (int)((options.zhigh - options.zlow) / options.dz);
    n = (simulation.n_points_.x + 1) * (simulation.n_points_.y + 1) * (simulation.n_points_.z + 1);
    debugline(simulation, simulation.root_processor_, "3D data grid is: %d x %d x %d (%d in total)",
              1 + simulation.n_points_.x, 1 + simulation.n_points_.y,
              1 + simulation.n_points_.z, n);

#if !SERIAL
    /* Wait! */
    MPI_Barrier(simulation.universe);
#endif

    /* If we which to precompute the force grid, do it now */
    // if (Options.rigidgrid) {
        // build3DForceGrid();
    // }

    /* Open all the file streams (one for every z point) [ONLY ON ROOT PROCESSOR] */
    if (simulation.onRootProcessor()) {
        simulation.fstreams_.reserve(simulation.n_points_.z + 1);
        for (int i = 0; i <= simulation.n_points_.z; ++i) {
            z = options.zhigh - i*options.dz;
            if (options.gzip == TRUE) {
                sprintf(outfile, "gzip -6 > scan-%06.3f.dat.gz", z);
                simulation.fstreams_.push_back(popen(outfile, "w"));
            }
            else {
                sprintf(outfile, "scan-%06.3f.dat", z);
                simulation.fstreams_.push_back(fopen(outfile, "w"));
            }
        }
    }

    /* Note the time */
    simulation.time_start_ = chrono::system_clock::now();
    return;
}

/* Now move the tip! */
// void moveTip(void) {

    // int n, nmax, check;
    // int i, ix, iy, iz;
    // double x, y, z;
    // double angle, minangle, maxforce, dd;
    // double e, ediff, eold, fnorm;
    // VECTOR f, d;
    // VECTOR *ftip;
    // int nxy, onproc, bufsize, nsr;
    // BUFFER *sendbuf, *recvbuf;
    // double checkperc, curperc;

    // FILE *tfp;
    // char dump[NAME_LENGTH];

    // [> Create a force vector (for real time analysis) <]
    // ftip = (VECTOR *)malloc((Npoints.z+1)*sizeof(VECTOR));

    // [> DEBUG DEBUG DEBUG DEBUG DEBUG <]
    // //Npoints.x = Npoints.y = 10;

    // [> Some initialization <]
    // Ntotal = 0;
    // checkperc = curperc = 0.10;
    // debugline(simulation, RootProc,"Simulation run started now");

    // [> Initialize the storage send and receive buffers <]
    // nxy = (Npoints.x + 1)*(Npoints.y + 1);
    // bufsize = Options.bufsize * (Npoints.z+1);
    // sendbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    // recvbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    // nsr = 0;

    // [> Loop x <]
    // for (ix=0; ix<=Npoints.x; ++ix) {
        // x = ix*Options.dx; [> Current x <]

        // [> Loop y <]
        // for (iy=0; iy<=Npoints.y; ++iy) {
            // y = iy*Options.dy; [> Current y <]

            // [> Check the progress and report every so often <]
            // n = ix*(Npoints.y+1) + iy;
            // if ( (Me == RootProc) && ((((double)n)/(nxy)) >= curperc) ) {
                // debugline(simulation, RootProc,"Finished approximately %4.1f %% of the simulation",100*curperc);
                // curperc += checkperc;
            // }

            // [> Compute on which processor this x,y combination should be run <]
            // onproc = n % NProcessors;
            // if (onproc != Me) {
                // continue;
            // }
            // PointsOnProc[Me]++;

            // [> Position the tip far above the surface <]
            // Dummy_pos.x = x;
            // Dummy_pos.y = y;
            // Dummy_pos.z = Options.zhigh + Options.dz; [> Plus dz to allow for initial subtraction <]
            // Tip_pos.x = Dummy_pos.x;
            // Tip_pos.y = Dummy_pos.y;
            // Tip_pos.z = Dummy_pos.z - DummyParams.rmin;

            // [> If the molecule is flexible, reset it to its original position, before beginning the approach <]
            // if (Options.flexible) {
                // for (i=0; i<Natoms; ++i) {
                    // Surf_pos[i].x = Surf_pos_org[i].x;
                    // Surf_pos[i].y = Surf_pos_org[i].y;
                    // Surf_pos[i].z = Surf_pos_org[i].z;
                // }
            // }

            // [> Approach and optimize <]
            // nmax = 0;
            // minangle = 9e99;
            // maxforce = -9e99;
            // for (iz=0; iz<=Npoints.z; ++iz) {
                // z = Options.zhigh - iz*Options.dz; [> Current z <]

                // [> Move tip and dummy atom toward the surface <]
                // Dummy_pos.z -= Options.dz;
                // Tip_pos.z -= Options.dz;

                // [> Collect the force <]
                // ftip[iz] = NULL_vector;

                // [> Relax/Minimize the configuration <]
                // ediff = 5*Options.etol;
                // for (n=0; n<Options.maxsteps; ++n) {
                    // nmax++;

                    // [> Compute all interaction energies <]
                    // interactTipSurface();
                    // interactTipDummy();
                    // interactTipHarmonic();

                    // [> Total energy and force computed <]
                    // e = TipSurf_energy + TipDummy_energy + TipHarmonic_energy;
                    // f.x = TipSurf_force.x + TipDummy_force.x + TipHarmonic_force.x;
                    // f.y = TipSurf_force.y + TipDummy_force.y + TipHarmonic_force.y;
                    // f.z = TipSurf_force.z + TipDummy_force.z + TipHarmonic_force.z;
                    // fnorm = sqrt(f.x*f.x + f.y*f.y + f.z*f.z);

                    // [> Energy difference <]
                    // if (n>0) {
                        // ediff = e - eold;
                    // }
                    // eold = e;

                    // [> Are the forces/energies tolerable <]
                    // if (Options.minterm == MIN_E) {
                        // check = (fabs(ediff)<Options.etol);
                    // }
                    // else if (Options.minterm == MIN_F) {
                        // check = (fabs(fnorm)<Options.ftol);
                    // }
                    // else if (Options.minterm == MIN_EF) {
                        // check = ((fabs(ediff)<Options.etol)&&(fabs(fnorm)<Options.ftol));
                    // }
                    // if (check) {
                        // break;
                    // }

                    // [> If they are not tolerable, update position of the tip atom based on the force <]
                    // Tip_pos.x += Options.cfac * f.x;
                    // Tip_pos.y += Options.cfac * f.y;
                    // Tip_pos.z += Options.cfac * f.z;

                    // [> If the molecule is supposed to behave flexible, we need to update its positions too <]
                    // if (Options.flexible) {
                        // updateFlexibleMolecule();
                    // }

                // } [> End minimization loop <]

                // [> Compute some other interesting data <]
                // d.x = Tip_pos.x - x;
                // d.y = Tip_pos.y - y;
                // d.z = Tip_pos.z - z;
                // dd = sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
                // angle = atan2(sqrt(d.x*d.x + d.y*d.y),d.z)*(180.0/PI);
                // if (angle<minangle) {
                    // minangle = angle;
                // }
                // if (fnorm>maxforce) {
                    // maxforce = fnorm;
                // }

                // [> Collect the final forces on the tip because of the surface <]
                // ftip[iz].x = TipSurf_force.x;
                // ftip[iz].y = TipSurf_force.y;
                // ftip[iz].z = TipSurf_force.z;

                // [> Store data in send buffers <]
                // sendbuf[nsr].ix = ix;
                // sendbuf[nsr].iy = iy;
                // sendbuf[nsr].iz = iz;
                // sendbuf[nsr].n  = n;
                // sendbuf[nsr].pos.x = x;
                // sendbuf[nsr].pos.y = y;
                // sendbuf[nsr].pos.z = z;
                // sendbuf[nsr].f.x = TipSurf_force.x;
                // sendbuf[nsr].f.y = TipSurf_force.y;
                // sendbuf[nsr].f.z = TipSurf_force.z;
                // sendbuf[nsr].d.x = d.x;
                // sendbuf[nsr].d.y = d.y;
                // sendbuf[nsr].d.z = d.z;
                // sendbuf[nsr].dd = dd;
                // sendbuf[nsr].e = TipSurf_energy;
                // sendbuf[nsr].angle = angle;
                // nsr++;

            // } [> End loop in z <]

            // [> Dump to file (only when the buffer is full, this can and should only happen after a full z approach) <]
            // if (nsr==bufsize) {
                // dumpToFiles(sendbuf,recvbuf,bufsize);
                // nsr = 0;
            // }

            // [> Compute the frequency shift for the given F(z) <]
            // //computeDeltaF(x,y,ftip);

            // if ((ix==(Npoints.x/2)) && (iy==(Npoints.y/2))) {
                // sprintf(dump,"test-%d-%d.xyz",ix,iy);
                // tfp = fopen(dump,"w");
                // fprintf(tfp,"%d\n\n",Natoms);
                // for (i=0; i<Natoms; ++i) {
                    // fprintf(tfp,"%s %8.4f %8.4f %8.4f\n",Surf_type[i],Surf_pos[i].x,Surf_pos[i].y,Surf_pos[i].z);
                // }
                // fclose(tfp);
            // }

            // [> Keep track of counting <]
            // Ntotal += nmax;

        // } [> End loop in y <]

    // } [> End loop in x <]

    // [> Dump to file (if it happened that the buffer contains anything at all <]
// #if !SERIAL
    // MPI_Allreduce(&nsr,&n,1,MPI_INT,MPI_SUM,Universe);
// #else
    // n = nsr;
// #endif
    // if (n>0) {
        // dumpToFiles(sendbuf,recvbuf,nsr);
    // }

    // [> Say one last thing <]
    // debugline(simulation, RootProc,"Finished %5.1f %% of the simulation",100*curperc);

    // [> Go home <]
    // return;
// }

/* Close all the file streams */
void closeUniverse(Simulation& simulation) {

#if !SERIAL
    /* Wait! */
    MPI_Barrier(simulation.universe);
#endif

    /* Close each separate file stream [ONLY ON ROOT PROCESSOR] */
    if (simulation.onRootProcessor()) {
        for (auto& file : simulation.fstreams_) {
            if (simulation.options_.gzip == TRUE) {
                pclose(file);
            }
            else {
                fclose(file);
            }
        }
    }
    return;
}

/********************
 ** FINAL THOUGHTS **
 ********************/

void finalize(Simulation& simulation) {

    int n, nsum;
    chrono::duration<double> dtime, timesum;

    /* Note the time */
    simulation.time_end_ = chrono::system_clock::now();

    /* Time difference */
    dtime = simulation.time_end_ - simulation.time_start_;
    timesum = timesum.zero();
#if !SERIAL
    MPI_Reduce(&dtime, &timesum, 1, MPI_DOUBLE, MPI_SUM, simulation.root_processor_,
                                                         simulation.universe);
#else
    timesum += dtime;
#endif

    /* Collect number of steps from all processors */
    nsum = 0;
#if !SERIAL
    MPI_Reduce(&simulation.n_total_, &nsum, 1, MPI_INT, MPI_SUM, simulation.root_processor_,
                                                                 simulation.universe);
#else
    nsum += simulation.n_total_;
#endif

    /* Print some miscelleneous information */
    debugline(simulation, simulation.root_processor_, "Simulation run finished");
    debugline(simulation, simulation.root_processor_, "Statistics:");
    n = (simulation.n_points_.x + 1) * (simulation.n_points_.y + 1) * (simulation.n_points_.z + 1);
    debugline(simulation, simulation.root_processor_, "    Computed %ld tip positions", n);
    debugline(simulation, simulation.root_processor_, "    Needed %ld minimization steps in total", nsum);
    debugline(simulation, simulation.root_processor_, "    Which means approximately %.2f minimization steps per tip position",
                                                                ((double) nsum / n));
    debugline(simulation, simulation.root_processor_, "    The simulation wall time is %.2f seconds", timesum);
    debugline(simulation, simulation.root_processor_, "    The entire simulation took %.2f seconds", dtime);
    debugline(simulation, simulation.root_processor_, "");
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
    MPI_Comm_rank(simulation.universe, &simulation.me_);
    MPI_Comm_size(simulation.universe, &simulation.n_processors_);
#else
    simulation.me_ = 0;
    simulation.n_processors_ = 1;
#endif
    /* Initialize the checker on how many x,y points for each processor */
    simulation.points_per_processor_.reserve(simulation.n_processors_);
    for (int i = 0; i < simulation.n_processors_; ++i) {
        simulation.points_per_processor_.push_back(0);
    }
    return;
}

/* Terminate our parallel worlds */
void closeParallelUniverse(Simulation& simulation) {

    vector<int> pop;

    /* How many x,y points on each processor */
#if !SERIAL
    MPI_Barrier(simulation.universe);
#endif
    pop.reserve(simulation.n_processors_);
    for (int i = 0; i < simulation.n_processors_; ++i) {
        pop.push_back(0);
    }
#if !SERIAL
    MPI_Reduce(static_cast<void*>(simulation.points_per_processor_.data()),
               static_cast<void*>(pop.data()), simulation.n_processors_,
               MPI_INT, MPI_SUM, simulation.root_processor_, simulation.universe);
#else
    pop[simulation.me_] += simulation.points_per_processor_[simulation.me_];
#endif
    debugline(simulation, simulation.root_processor_, "How many x,y points did each processor handle:");
    for (int i = 0; i < simulation.n_processors_; ++i) {
        debugline(simulation, simulation.root_processor_, "    Processor %2d: %6d x,y points", i, pop[i]);
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
    parseCommandLine(argc, argv, simulation);
    readInputFile(simulation);
    readXYZFile(simulation);
    readParameterFile(simulation);

    /* If the molecule is flexible, build the topology */
    // if (Options.flexible) {
        // buildTopology();
    // }

    /* The simulation itself */
    openUniverse(simulation);
    //moveTip();
    closeUniverse(simulation);

    /* Some final thoughts */
    finalize(simulation);

    /* And stop the parallel routines properly */
    closeParallelUniverse(simulation);
    return 0;
}
