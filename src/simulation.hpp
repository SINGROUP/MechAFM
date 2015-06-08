#pragma once

#include <vector>

#include "system.hpp"
#include "vectors.hpp"

using namespace std;

/* Define a structure for all input options */
typedef struct InputOptions {
    char xyzfile[NAME_LENGTH];
    char paramfile[NAME_LENGTH];
    char tipatom[NAME_LENGTH];
    char dummyatom[NAME_LENGTH];
    char planeatom[NAME_LENGTH];
    double dx, dy, dz;
    double zlow, zhigh, zplane;
    int coulomb, units;
    int maxsteps, minterm;
    double etol, ftol, cfac;
    int bufsize, gzip;
    int flexible, rigidgrid;
} InputOptions;

class Simulation {
 public:
    Simulation(): {};
    ~Simulation() {};
    System system;
    char input_file_name_[NAME_LENGTH];        /* File name of the input file */
    InputOptions options_;                   /* Structure containing all relevant input options */
    Vec3d box_;                             /* Vector containing the size of the universe */
    char **SurfType2Num;                    /* Dictionary hash to go from atom types to numbers in the 2D array (flexible molecule) */
    FILE **FStreams;                        /* Array with the entire file stream */

    /* Some grid computing thingies */
    double *ForceGridRigid;                 /* 3D force grid for use with rigid tips (interpolation) */
    Vec3i Ngrid;                          /* Size of 3D force grid */
    int Ngridpoints;                        /* Total number of gridpoints */
    double GridSpacing;                     /* The size of the cubes of the grid */

    /* Some parallel specific global variables */
#if !SERIAL
        MPI_Comm universe;                  /* The entire parallel universe */
#endif
    int n_processors_;                        /* Total number of processors */
    int me_;                                 /* The current processor */
    int root_processor_;                           /* The main processor */
    vector<int> points_per_processor;       /* How many x,y points on this processor */

struct timeval TimeStart, TimeEnd;


 private:
};
