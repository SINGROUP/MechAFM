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
#if MPI_BUILD
    #include <mpi.h>
#endif
#include <chrono>

#include "globals.hpp"
#include "messages.hpp"
#include "parse.hpp"
#include "simulation.hpp"
#include "system.hpp"
#include "interactions.hpp"

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
    simulation.n_total_ = 0;
    simulation.n_points_.x = (int)(simulation.box_.x / options.dx);
    simulation.n_points_.y = (int)(simulation.box_.y / options.dy);
    simulation.n_points_.z = (int)((options.zhigh - options.zlow) / options.dz);
    n = (simulation.n_points_.x + 1) * (simulation.n_points_.y + 1) * (simulation.n_points_.z + 1);
    pretty_print(simulation, "3D data grid is: %d x %d x %d (%d in total)",
              1 + simulation.n_points_.x, 1 + simulation.n_points_.y,
              1 + simulation.n_points_.z, n);

#if MPI_BUILD
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

/* Close all the file streams */
void closeUniverse(Simulation& simulation) {

#if MPI_BUILD
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
#if MPI_BUILD
    MPI_Reduce(&dtime, &timesum, 1, MPI_DOUBLE, MPI_SUM, simulation.root_processor_,
                                                         simulation.universe);
#else
    timesum += dtime;
#endif

    /* Collect number of steps from all processors */
    nsum = 0;
#if MPI_BUILD
    MPI_Reduce(&simulation.n_total_, &nsum, 1, MPI_INT, MPI_SUM, simulation.root_processor_,
                                                                 simulation.universe);
#else
    nsum += simulation.n_total_;
#endif

    /* Print some miscelleneous information */
    pretty_print(simulation, "Simulation run finished");
    pretty_print(simulation, "Statistics:");
    n = (simulation.n_points_.x + 1) * (simulation.n_points_.y + 1) * (simulation.n_points_.z + 1);
    pretty_print(simulation, "    Computed %ld tip positions", n);
    pretty_print(simulation, "    Needed %ld minimization steps in total", nsum);
    pretty_print(simulation, "    Which means approximately %.2f minimization steps per tip position",
                                                                ((double) nsum / n));
    pretty_print(simulation, "    The simulation wall time is %.2f seconds", timesum);
    pretty_print(simulation, "    The entire simulation took %.2f seconds", dtime);
    pretty_print(simulation, "");
    return;
}

/***************************
 ** THE PARALLEL UNIVERSE **
 ***************************/

/* Initialize our parallel world */
void openParallelUniverse(int argc, char *argv[], Simulation& simulation) {
#if MPI_BUILD
    /* Start MPI */
    MPI_Init(&argc,&argv);
#endif
    /* Determine the size of the universe and which processor we are on */
    simulation.root_processor_ = 0;
#if MPI_BUILD
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
#if MPI_BUILD
    MPI_Barrier(simulation.universe);
#endif
    pop.reserve(simulation.n_processors_);
    for (int i = 0; i < simulation.n_processors_; ++i) {
        pop.push_back(0);
    }
#if MPI_BUILD
    MPI_Reduce(static_cast<void*>(simulation.points_per_processor_.data()),
               static_cast<void*>(pop.data()), simulation.n_processors_,
               MPI_INT, MPI_SUM, simulation.root_processor_, simulation.universe);
#else
    pop[simulation.me_] += simulation.points_per_processor_[simulation.me_];
#endif
    pretty_print(simulation, "How many x,y points did each processor handle:");
    for (int i = 0; i < simulation.n_processors_; ++i) {
        pretty_print(simulation, "    Processor %2d: %6d x,y points", i, pop[i]);
    }

#if MPI_BUILD
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

    simulation.buildInteractions();

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
