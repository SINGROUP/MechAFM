#pragma once

#include <string>
#include <vector>

#include "globals.hpp"
#include "vectors.hpp"

using namespace std;

/* Define a structure to quickly compute nonbonded interactions */
typedef struct InteractionList {
    double eps;
    double sig;
    double es12;
    double es6;
    double qq;
    double k;
    double r0;
    double rmin;
    int morse;
    double De;
    double a;
    double re;
} InteractionList;

/* Define a structure to figure out possible bonds */
typedef struct PosBonds {
    char a1[NAME_LENGTH], a2[NAME_LENGTH];
    double r0;
} PosBonds;

/* Define a structure to compute bonds */
typedef struct BondInteraction {
    int a1, a2;
    double r0;
    double k;
} BondInteraction;

/* Define a structure to compute angles */
typedef struct AngleInteraction {
    int a1, a2, a3;
    int bond1, bond2;
    double theta0;
    double k;
} AngleInteraction;

class System {
 public:
     System() {};
     ~System() {};
    void initialize(int n_atoms);
    inline void setTipDummyDistance(double d) {tip_dummy_d_ = d;}
    void setDummyXY(double x, double y);
    void setDummyZ(double z);

    int n_atoms_;                             /* Number of surface atoms */
    int n_types_;                             /* Number of surface atom types (flexible molecule) */

    // Vectors holding the system state
    // index 0 = dummy and index 1 = tip
    vector<Vec3d> positions_;
    vector<Vec3d> velocities_;
    vector<Vec3d> forces_;
    vector<double> charges_;
    vector<double> masses_;
    vector<int> fixed_;
    vector<string> types_;

    Vec3d TipSurf_force;                   /* Force on tip atom caused by the surface */
    Vec3d TipDummy_force;                  /* Force on tip atom caused by the dummy atom */
    Vec3d TipHarmonic_force;               /* Force on tip atom caused by the harmonic constraint */
    double TipSurf_energy;                  /* Energy of tip atom caused by the surface */
    double TipDummy_energy;                 /* Energy of tip atom caused by the dummy atom */
    double TipHarmonic_energy;              /* Energy of tip atom caused by the harmonic constraint */
    // VECTOR *Surf_pos;                       [> Vector containing positions of all surface atoms <]
    // VECTOR *Surf_vel;                       [> Vector containing velocities of all surface atoms (needed for flexible/FIRE minimizer) <]
    // VECTOR *Surf_pos_org;                   [> Copy of the surface position vector, needed for flexible reset <]
    // VECTOR *Surf_force;                     [> Vector containing the surface forces (needed for flexible/FIRE minimizer) <]
    // double *Surf_q;                         [> List containing charges of all surface atoms <]
    // double *Surf_mass;                      [> List containing masses of all surface atoms <]
    // int *Surf_fix;                          [> List containing a boolean to signal fixed or free atoms (needed for flexible) <]
    // char **Surf_type;                       [> List containing types of all surface atoms <]
    // InteractionList *TipSurfParams;         [> Structured list for all tip-surface particle interaction parameters <]
    // InteractionList DummyParams;            [> List for particle interaction parameters with dummy atom <]
    // InteractionList Harmonic;               [> List for the harmonic constraint parameters on the tip atom <]
    // InteractionList **SurfSurfParams;       [> 2D array for all surface-surface particle interaction parameters (flexible molecule) <]
    // BondInteraction *Bonds;                 [> List of all possible bonds (flexible molecule) <]
    // AngleInteraction *Angles;               [> List of all possible angles (flexible molecule) <]
    // InteractionList Substrate;              [> Substrate support parameters (flexible molecule) <]

 private:
    double tip_dummy_d_;
};
