#pragma once

#if MPI_BUILD
    #include <mpi.h>
#endif
#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include "globals.hpp"
#include "integrators.hpp"
#include "interactions.hpp"
#include "minimiser.hpp"
#include "system.hpp"
#include "vectors.hpp"

using namespace std;

/* A simple list to distinguish different minimization criteria */
enum MinimizationCriteria {MIN_E, MIN_F, MIN_EF, NOT_SET};

/* A simple list to distinguish the chosen unit system */
enum Units {U_KCAL, U_KJ, U_EV};

/* Define a structure for all input options */
struct InputOptions {
    char xyzfile[NAME_LENGTH];
    char paramfile[NAME_LENGTH];
    char tipatom[NAME_LENGTH];
    char dummyatom[NAME_LENGTH];
    char planeatom[NAME_LENGTH];
    Vec3d box;                     // Vector containing the size of the universe
    double dx, dy, dz;
    double zlow, zhigh, zplane;
    Units units;
    bool coulomb;
    int maxsteps;
    MinimizationCriteria minterm;
    double etol, ftol, cfac;
    int bufsize;
    bool gzip;
    bool flexible, rigidgrid;
    bool xyz_charges;
    MinimiserType minimiser_type;
    IntegratorType integrator_type;
};

class Simulation {
 public:
    Simulation() {};
    ~Simulation() {};
    bool rootProcess();
    void run();
    void buildInteractions();

    System system;
    char input_file_name_[NAME_LENGTH];     /* File name of the input file */
    InputOptions options_;                  /* Structure containing all relevant input options */
    InteractionParameters interaction_parameters_;
    Vec3i n_points_;                        /* Number of points (x,y,z) for the tip */
    long int n_total_;                      /* Total number of minimization loops used */
    vector<FILE*> fstreams_;                /* Array with the entire file stream */
    vector<unique_ptr<Interaction>> interactions_;

    /* Some grid computing thingies */
    double *ForceGridRigid;                 /* 3D force grid for use with rigid tips (interpolation) */
    Vec3i Ngrid;                            /* Size of 3D force grid */
    int Ngridpoints;                        /* Total number of gridpoints */
    double GridSpacing;                     /* The size of the cubes of the grid */

    /* Some parallel specific global variables */
#if MPI_BUILD
        MPI_Comm universe;                  /* The entire parallel universe */
#endif
    int n_processes_;                      /* Total number of processes */
    int root_process_;
    int current_process_;                                /* The current process */
    vector<int> points_per_process_;      /* How many x,y points on this process */

    chrono::time_point<chrono::system_clock> time_start_, time_end_;

 private:
    void calculateTipDummyDistance();
    void writeOutput(vector<OutputData> output_buffer);
    void addLJInteraction(int atom_i1, int atom_i2);
    void addCoulombInteraction(int atom_i1, int atom_i2);
    bool findOverwriteParameters(int atom_i1, int atom_i2, OverwriteParameters op);
    void buildTipSurfaceInteractions();
    void buildTipGridInteractions();
    void buildTipDummyInteractions();
    void buildSurfaceSurfaceInteractions();
    void buildSubstrateInteractions();
};
