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

// Defines all the different minimization criteria
enum MinimizationCriteria {MIN_E, MIN_F, MIN_EF, NOT_SET};

// Defines the possible unit systems
enum Units {U_KCAL, U_KJ, U_EV};

// Defines a structure for all input options
struct InputOptions {
    string inputfolder;
    string outputfolder;
    string inputfile;
    string xyzfile;
    string paramfile;
    string tipatom;
    string dummyatom;
    string planeatom;
    Vec3d box;
    double dx, dy, dz;
    double zlow, zhigh, zplane;
    Units units;
    bool coulomb;
    int maxsteps;
    MinimizationCriteria minterm;
    double etol, ftol, dt;
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
    // Returns whether we're on the root process or not
    bool rootProcess();
    // Runs the simulation
    void run();
    // Builds all the interactions
    void buildInteractions();

    System system;  // Holds the system to be minimised
    vector<unique_ptr<Interaction>> interactions_; // List of all the interactions
    InputOptions options_;  // Structure containing all relevant input options
    InteractionParameters interaction_parameters_;
    Vec3i n_points_;  // Number of points (x,y,z) to be minimised
    long int n_total_;  // Total number of minimization steps used
    vector<FILE*> fstreams_;  // Array with all the file streams

    // Some parallel specific global variables
#if MPI_BUILD
        MPI_Comm universe;  //The entire parallel universe
#endif
    int n_processes_;  // Total number of processes
    int root_process_;  // The number of the root process (usually 0)
    int current_process_;  // The number of the current process
    vector<int> points_per_process_;  // How many x,y points on each process

    // Simulation start and end time
    chrono::time_point<chrono::system_clock> time_start_, time_end_;

 private:
    // Calculates the initial distance of the tip and the dummy atoms
    void calculateTipDummyDistance();
    // Writes the output buffer to the disk
    void writeOutput(vector<OutputData> output_buffer);
    // Add a LJ or Morse interaction between atoms 1 and 2
    void addVDWInteraction(int atom_i1, int atom_i2);
    // Add a Coulomb interaction between atoms 1 and 2
    void addCoulombInteraction(int atom_i1, int atom_i2);
    // Looks for overwrite parameters for atoms 1 and 2. Return true if found
    // and sets op if found.
    bool findOverwriteParameters(int atom_i1, int atom_i2, OverwriteParameters& op);
    // Build all the interactions of the tip atom with the surface atoms
    void buildTipSurfaceInteractions();
    // Build a grid interaction to approximate tip surface interactions
    void buildTipGridInteractions();
    // Build the interactions between the tip and dummy atom
    void buildTipDummyInteractions();
    // Build all the interactions between the surface atoms
    void buildSurfaceSurfaceInteractions();
    // Build substrate interactions for all surface atoms
    void buildSubstrateInteractions();
};
