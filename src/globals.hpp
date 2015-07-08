#pragma once

#include "vectors.hpp"

/* Some macro definitions */
#define LINE_LENGTH 512
#define NAME_LENGTH 128
#define ATOM_LENGTH 8
#define TOLERANCE 1e-10
#define TRUE 1
#define FALSE 0
#define NEGVAL -99999
#define PI 3.14159265358979323846
#define SIXTHRT2 1.12246204830937298142

/* Define a structure for easy processor communication */
struct OutputData {
    Vec3i indices;
    int minimisation_steps;
    double angle, tip_energy, r;
    Vec3d position, tip_force, r_vec;
};
