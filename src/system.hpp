#pragma once

#include <vector>

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
typedef struct PossibleBonds {
    char a1[NAME_LENGTH], a2[NAME_LENGTH];
    double r0;
} PossibleBonds;

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
     System(): {};
     ~System() {};
    int n_atoms_;                             /* Number of surface atoms */
    int Nbonds;                             /* Number of bonds in flexible molecule */
    int Nangles;                            /* Number of angles in flexible molecule */
    int n_types_;                             /* Number of surface atom types (flexible molecule) */
    int n_fixed_;                             /* Number of fixed atoms (flexible molecule) */
    // VECTOR *Surf_pos;                       [> Vector containing positions of all surface atoms <]
    // VECTOR *Surf_vel;                       [> Vector containing velocities of all surface atoms (needed for flexible/FIRE minimizer) <]
    // VECTOR *Surf_pos_org;                   [> Copy of the surface position vector, needed for flexible reset <]
    // VECTOR *Surf_force;                     [> Vector containing the surface forces (needed for flexible/FIRE minimizer) <]
    // double *Surf_q;                         [> List containing charges of all surface atoms <]
    // double *Surf_mass;                      [> List containing masses of all surface atoms <]
    // int *Surf_fix;                          [> List containing a boolean to signal fixed or free atoms (needed for flexible) <]
    // char **Surf_type;                       [> List containing types of all surface atoms <]
    vector<Vec3d> positions_;
    vector<Vec3d> velocities_;
    vector<Vec3d> forces_;
    vector<double> charges_;
    vector<double> masses_;
    vector<bool> fixed_;
    vector<string> types_;
    InteractionList *TipSurfParams;         /* Structured list for all tip-surface particle interaction parameters */
    InteractionList DummyParams;            /* List for particle interaction parameters with dummy atom */
    InteractionList Harmonic;               /* List for the harmonic constraint parameters on the tip atom */
    InteractionList **SurfSurfParams;       /* 2D array for all surface-surface particle interaction parameters (flexible molecule) */
    BondInteraction *Bonds;                 /* List of all possible bonds (flexible molecule) */
    AngleInteraction *Angles;               /* List of all possible angles (flexible molecule) */
    InteractionList Substrate;              /* Substrate support parameters (flexible molecule) */
    VECTOR Tip_pos;                         /* Position of the tip atom */
    VECTOR Tip_vel;                         /* Velocities of the tip atom */
    double Tip_mass;                        /* Mass of the tip atom */
    VECTOR Dummy_pos;                       /* Position of the dummy atom */
    VECTOR TipSurf_force;                   /* Force on tip atom caused by the surface */
    VECTOR TipDummy_force;                  /* Force on tip atom caused by the dummy atom */
    VECTOR TipHarmonic_force;               /* Force on tip atom caused by the harmonic constraint */
    double TipSurf_energy;                  /* Energy of tip atom caused by the surface */
    double TipDummy_energy;                 /* Energy of tip atom caused by the dummy atom */
    double TipHarmonic_energy;              /* Energy of tip atom caused by the harmonic constraint */
    Vec3i Npoints;                        /* Number of points (x,y,z) for the tip */
    long int Ntotal;                        /* Total number of minimization loops used */

 private:
};
