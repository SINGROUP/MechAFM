#pragma once

#if MPI_BUILD
    #include <mpi.h>
#endif
#include <chrono>
#include <string>
#include <vector>

#include "globals.hpp"
#include "interactions.hpp"
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
    double dx, dy, dz;
    double zlow, zhigh, zplane;
    Units units;
    bool coulomb;
    int maxsteps;
    MinimizationCriteria minterm;
    double etol, ftol, cfac;
    int bufsize, gzip;
    bool flexible, rigidgrid;
    bool xyz_charges;
};

class Simulation {
 public:
    Simulation() {};
    ~Simulation() {};
    inline bool onRootProcessor() {return me_ == root_processor_;}

    System system;
    char input_file_name_[NAME_LENGTH];        /* File name of the input file */
    InputOptions options_;                   /* Structure containing all relevant input options */
    InteractionParameters interaction_parameters_;
    Vec3d box_;                             /* Vector containing the size of the universe */
    Vec3i n_points_;                        /* Number of points (x,y,z) for the tip */
    long int n_total_;                        /* Total number of minimization loops used */
    char **SurfType2Num;                    /* Dictionary hash to go from atom types to numbers in the 2D array (flexible molecule) */
    vector<FILE*> fstreams_;                        /* Array with the entire file stream */

    /* Some grid computing thingies */
    double *ForceGridRigid;                 /* 3D force grid for use with rigid tips (interpolation) */
    Vec3i Ngrid;                          /* Size of 3D force grid */
    int Ngridpoints;                        /* Total number of gridpoints */
    double GridSpacing;                     /* The size of the cubes of the grid */

    /* Some parallel specific global variables */
#if MPI_BUILD
        MPI_Comm universe;                  /* The entire parallel universe */
#endif
    int n_processors_;                        /* Total number of processors */
    int me_;                                 /* The current processor */
    int root_processor_;                           /* The main processor */
    vector<int> points_per_processor_;       /* How many x,y points on this processor */

    chrono::time_point<chrono::system_clock> time_start_, time_end_;

 private:
};
