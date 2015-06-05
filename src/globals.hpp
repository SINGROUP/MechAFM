#pragma once

#if !SERIAL
    #include <mpi.h>
#endif

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

/* Local definitions */

/* Define vector types */
typedef struct vector {
    double x, y, z;
} VECTOR;

typedef struct ivector {
    int x, y, z;
} IVECTOR;

extern IVECTOR NULL_ivector;
extern VECTOR NULL_vector;


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

/* Define a structure for easy processor communication */
typedef struct buffer {
    int ix, iy, iz, n;
    double angle, e, dd;
    VECTOR pos, f, d;
} BUFFER;

/* A simple list to distinguish different minimization criteria */
enum{MIN_E,MIN_F,MIN_EF};

/* A simple list to distinguish the chosen unit system */
enum{U_KCAL,U_KJ,U_EV};

/* Time thingies */
extern struct timeval TimeStart, TimeEnd;

/**********************
 ** GLOBAL VARIABLES **
 **********************/

extern char InputFileName[NAME_LENGTH];        /* File name of the input file */
extern InputOptions Options;                   /* Structure containing all relevant input options */
extern int Natoms;                             /* Number of surface atoms */
extern int Nbonds;                             /* Number of bonds in flexible molecule */
extern int Nangles;                            /* Number of angles in flexible molecule */
extern int Ntypes;                             /* Number of surface atom types (flexible molecule) */
extern int Nfixed;                             /* Number of fixed atoms (flexible molecule) */
extern VECTOR *Surf_pos;                       /* Vector containing positions of all surface atoms */
extern VECTOR *Surf_vel;                       /* Vector containing velocities of all surface atoms (needed for flexible/FIRE minimizer) */
extern VECTOR *Surf_pos_org;                   /* Copy of the surface position vector, needed for flexible reset */
extern VECTOR *Surf_force;                     /* Vector containing the surface forces (needed for flexible/FIRE minimizer) */
extern double *Surf_q;                         /* List containing charges of all surface atoms */
extern double *Surf_mass;                      /* List containing masses of all surface atoms */
extern int *Surf_fix;                          /* List containing a boolean to signal fixed or free atoms (needed for flexible) */
extern char **Surf_type;                       /* List containing types of all surface atoms */
extern VECTOR Box;                             /* Vector containing the size of the universe */
extern InteractionList *TipSurfParams;         /* Structured list for all tip-surface particle interaction parameters */
extern InteractionList DummyParams;            /* List for particle interaction parameters with dummy atom */
extern InteractionList Harmonic;               /* List for the harmonic constraint parameters on the tip atom */
extern InteractionList **SurfSurfParams;       /* 2D array for all surface-surface particle interaction parameters (flexible molecule) */
extern char **SurfType2Num;                    /* Dictionary hash to go from atom types to numbers in the 2D array (flexible molecule) */
extern BondInteraction *Bonds;                 /* List of all possible bonds (flexible molecule) */
extern AngleInteraction *Angles;               /* List of all possible angles (flexible molecule) */
extern InteractionList Substrate;              /* Substrate support parameters (flexible molecule) */
extern VECTOR Tip_pos;                         /* Position of the tip atom */
extern VECTOR Tip_vel;                         /* Velocities of the tip atom */
extern double Tip_mass;                        /* Mass of the tip atom */
extern VECTOR Dummy_pos;                       /* Position of the dummy atom */
extern VECTOR TipSurf_force;                   /* Force on tip atom caused by the surface */
extern VECTOR TipDummy_force;                  /* Force on tip atom caused by the dummy atom */
extern VECTOR TipHarmonic_force;               /* Force on tip atom caused by the harmonic constraint */
extern double TipSurf_energy;                  /* Energy of tip atom caused by the surface */
extern double TipDummy_energy;                 /* Energy of tip atom caused by the dummy atom */
extern double TipHarmonic_energy;              /* Energy of tip atom caused by the harmonic constraint */
extern IVECTOR Npoints;                        /* Number of points (x,y,z) for the tip */
extern long int Ntotal;                        /* Total number of minimization loops used */
extern FILE **FStreams;                        /* Array with the entire file stream */

/* Some grid computing thingies */
extern double *ForceGridRigid;                 /* 3D force grid for use with rigid tips (interpolation) */
extern IVECTOR Ngrid;                          /* Size of 3D force grid */
extern int Ngridpoints;                        /* Total number of gridpoints */
extern double GridSpacing;                     /* The size of the cubes of the grid */

/* Function pointers */
extern void (*interactTipSurface)(void);       // For the tip surface interaction
extern void interactTipSurfaceDirectly(void);
extern void interactTipSurfaceFromGrid(void);

/* Some parallel specific global variables */
#if !SERIAL
    extern MPI_Comm Universe;                  /* The entire parallel universe */
#endif
extern int NProcessors;                        /* Total number of processors */
extern int Me;                                 /* The current processor */
extern int RootProc;                           /* The main processor */
extern int *PointsOnProc;                      /* How many x,y points on this processor */
