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
    simulation.n_points_.x = floor(simulation.options_.box.x / options.dx) + 1;
    simulation.n_points_.y = floor(simulation.options_.box.y / options.dy) + 1;
    simulation.n_points_.z = floor((options.zhigh - options.zlow) / options.dz) + 1;
    n = simulation.n_points_.x * simulation.n_points_.y * simulation.n_points_.z;
    pretty_print("3D data grid is: %d x %d x %d (%d in total)",
              simulation.n_points_.x, simulation.n_points_.y,
              simulation.n_points_.z, n);

#if MPI_BUILD
    /* Wait! */
    MPI_Barrier(simulation.universe);
#endif

    /* Open all the file streams (one for every z point) [ONLY ON ROOT PROCESS] */
    if (simulation.rootProcess()) {
        simulation.fstreams_.reserve(simulation.n_points_.z + 1);
        for (int i = 0; i <= simulation.n_points_.z - 1; ++i) {
            z = options.zhigh - i*options.dz;
            if (options.gzip == TRUE) {
                sprintf(outfile, "examples/gzip -6 > scan-%06.3f.dat.gz", z);
                simulation.fstreams_.push_back(popen(outfile, "w"));
            }
            else {
                sprintf(outfile, "examples/scan-%06.3f.dat", z);
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

    /* Close each separate file stream [ONLY ON ROOT PROCESS] */
    if (simulation.rootProcess()) {
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
    MPI_Reduce(&dtime, &timesum, 1, MPI_DOUBLE, MPI_SUM, simulation.root_process_,
                                                         simulation.universe);
#else
    timesum += dtime;
#endif

    /* Collect number of steps from all processes */
    nsum = 0;
#if MPI_BUILD
    MPI_Reduce(&simulation.n_total_, &nsum, 1, MPI_INT, MPI_SUM, simulation.root_process_,
                                                                 simulation.universe);
#else
    nsum += simulation.n_total_;
#endif

    /* Print some miscelleneous information */
    pretty_print("Simulation run finished");
    pretty_print("Statistics:");
    n = (simulation.n_points_.x) * (simulation.n_points_.y) * (simulation.n_points_.z);
    pretty_print("    Computed %ld tip positions", n);
    pretty_print("    Needed %ld minimization steps in total", nsum);
    pretty_print("    Which means approximately %.2f minimization steps per tip position",
                                                                ((double) nsum / n));
    pretty_print("    The simulation wall time is %.2f seconds", timesum);
    pretty_print("    The entire simulation took %.2f seconds", dtime);
    pretty_print("");
    return;
}

/***************************
 ** THE PARALLEL UNIVERSE **
 ***************************/

/* Initialize our parallel world */
void openParallelUniverse(int argc, char *argv[], Simulation& simulation) {
    (void)argc;
    (void)argv;
#if MPI_BUILD
    /* Start MPI */
    MPI_Init(&argc,&argv);
#endif
    /* Determine the size of the universe and which process we are on */
    simulation.root_process_ = 0;
#if MPI_BUILD
    simulation.universe = MPI_COMM_WORLD;
    MPI_Comm_rank(simulation.universe, &simulation.current_process_);
    MPI_Comm_size(simulation.universe, &simulation.n_processes_);
#else
    simulation.current_process_ = 0;
    simulation.n_processes_ = 1;
#endif
    /* Initialize the checker on how many x,y points for each process */
    simulation.points_per_process_.reserve(simulation.n_processes_);
    for (int i = 0; i < simulation.n_processes_; ++i) {
        simulation.points_per_process_.push_back(0);
    }
    return;
}

/* Terminate our parallel worlds */
void closeParallelUniverse(Simulation& simulation) {

    /* How many x,y points on each process */
#if MPI_BUILD
    MPI_Barrier(simulation.universe);
#endif
    vector<int> pop(simulation.n_processes_);
#if MPI_BUILD
    MPI_Reduce(static_cast<void*>(simulation.points_per_process_.data()),
               static_cast<void*>(pop.data()), simulation.n_processes_,
               MPI_INT, MPI_SUM, simulation.root_process_, simulation.universe);
#else
    pop[simulation.current_process_] += simulation.points_per_process_[simulation.current_process_];
#endif
    pretty_print("How many x,y points did each process handle:");
    for (int i = 0; i < simulation.n_processes_; ++i) {
        pretty_print("    Process %2d: %6d x,y points", i, pop[i]);
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

    /* The simulation itself */
    openUniverse(simulation);
    simulation.buildInteractions();
    simulation.run();
    closeUniverse(simulation);

    /* Some final thoughts */
    finalize(simulation);

    /* And stop the parallel routines properly */
    closeParallelUniverse(simulation);
    return 0;
}
