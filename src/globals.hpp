#pragma once

#define DEBUG_MODE false

#include "vectors.hpp"

// Some constant definitions
#define LINE_LENGTH 512
#define NAME_LENGTH 128
#define ATOM_LENGTH 8
#define TOLERANCE 1e-10
#define NEGVAL -99999
#define PI 3.14159265358979323846
#define SIXTHRT2 1.12246204830937298142

const double g_force_grid_margin = 1.5; // How wide of a margin force grid has around the simulation area when 'rigidgrid' is used
const double g_gaussian_cutoff_value = 1.0e-10; // Relative value of a Gaussian after which the rest of the values further away are approximated to zero
const double g_tip_gaussian_width = 0.5; // Width of the Gaussian charge distribution at the tip, in Å

// unit conversion factors
const double g_hartree_to_eV = 27.211386;
const double g_hartree_to_kJ = 2625.5002;
const double g_hartree_to_kcal = 627.50961;


// Define a structure for easy processor communication
struct OutputData {
    Vec3i indices;
    int minimisation_steps;
    double angle, tip_energy, r;
    Vec3d position, tip_force, r_vec;
};
