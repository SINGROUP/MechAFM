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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <glob.h>
#include <sys/time.h>

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

VECTOR NULL_vector = {0.0, 0.0, 0.0};

typedef struct ivector {
    int x, y, z;
} IVECTOR;

IVECTOR NULL_ivector = {0, 0, 0};

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
struct timeval TimeStart, TimeEnd;

/**********************
 ** GLOBAL VARIABLES **
 **********************/

char InputFileName[NAME_LENGTH];        /* File name of the input file */
InputOptions Options;                   /* Structure containing all relevant input options */
int Natoms;                             /* Number of surface atoms */
int Nbonds;                             /* Number of bonds in flexible molecule */
int Nangles;                            /* Number of angles in flexible molecule */
int Ntypes;                             /* Number of surface atom types (flexible molecule) */
int Nfixed;                             /* Number of fixed atoms (flexible molecule) */
VECTOR *Surf_pos;                       /* Vector containing positions of all surface atoms */
VECTOR *Surf_vel;                       /* Vector containing velocities of all surface atoms (needed for flexible/FIRE minimizer) */
VECTOR *Surf_pos_org;                   /* Copy of the surface position vector, needed for flexible reset */
VECTOR *Surf_force;                     /* Vector containing the surface forces (needed for flexible/FIRE minimizer) */
double *Surf_q;                         /* List containing charges of all surface atoms */
double *Surf_mass;                      /* List containing masses of all surface atoms */
int *Surf_fix;                          /* List containing a boolean to signal fixed or free atoms (needed for flexible) */
char **Surf_type;                       /* List containing types of all surface atoms */
VECTOR Box;                             /* Vector containing the size of the universe */
InteractionList *TipSurfParams;         /* Structured list for all tip-surface particle interaction parameters */
InteractionList DummyParams;            /* List for particle interaction parameters with dummy atom */
InteractionList Harmonic;               /* List for the harmonic constraint parameters on the tip atom */
InteractionList **SurfSurfParams;       /* 2D array for all surface-surface particle interaction parameters (flexible molecule) */
char **SurfType2Num;                    /* Dictionary hash to go from atom types to numbers in the 2D array (flexible molecule) */
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
IVECTOR Npoints;                        /* Number of points (x,y,z) for the tip */
long int Ntotal;                        /* Total number of minimization loops used */
FILE **FStreams;                        /* Array with the entire file stream */

/* Some grid computing thingies */
double *ForceGridRigid;                 /* 3D force grid for use with rigid tips (interpolation) */
IVECTOR Ngrid;                          /* Size of 3D force grid */
int Ngridpoints;                        /* Total number of gridpoints */
double GridSpacing;                     /* The size of the cubes of the grid */

/* Function pointers */
void (*interactTipSurface)(void);       /* For the tip surface interaction */
void interactTipSurfaceDirectly(void);
void interactTipSurfaceFromGrid(void);

/* Some parallel specific global variables */
#if !SERIAL
    MPI_Comm Universe;                  /* The entire parallel universe */
#endif
int NProcessors;                        /* Total number of processors */
int Me;                                 /* The current processor */
int RootProc;                           /* The main processor */
int *PointsOnProc;                      /* How many x,y points on this processor */

/**********************
 ** HEADER FUNCTIONS **
 **********************/

/* An error function */
void error(char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    fprintf(stderr, "+- ERROR (on proc %d): %s\n", Me, ws);
#if !SERIAL
    MPI_Finalize();
#endif
    exit(1);
}

/* And a warning function */
void warning(char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    if (Me==RootProc) {
        fprintf(stderr, "+- WARNING: %s\n", ws);
    }
}

/* Print some debug information to the screen */
void debugline(int proc, char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    if (Me==proc) {
        fprintf(stdout, "+- %s\n", ws);
    }
}

/* Convert a string to uppercase */
char *strupp(char *string) {
    char *convert;
    convert = string;
    do {
        *convert = toupper((unsigned char)*convert);
    } while (*convert++);
    return string;
}

/* Convert a string to lowercase */
char *strlow(char *string) {
    char *convert;
    convert = string;
    do {
        *convert = tolower((unsigned char)*convert);
    } while (*convert++);
    return string;
}

/* Filter out comment and empty lines in a file */
int checkForComments(char *line) {
    int moveon = FALSE;
    if (line[0] == '#') {
        moveon = TRUE;
    }
    else if (line[0] == '%') {
        moveon = TRUE;
    }
    else if (line[0] == '\n') {
        moveon = TRUE;
    }
    return moveon;
}

/* Check if a value is an integer */
int isint(char *str) {
    int integer = TRUE;
    int n = strlen(str);
    int i;
    for (i=0; i<n; ++i) {
        if (isdigit(str[i]) || str[i] == '-' || str[i] == '+'){
            continue;
        }
        integer = FALSE;
    }
    return integer;
}

/**************************
 ** FILE INPUT FUNCTIONS **
 **************************/

/* Read stuff from the command line */
void parseCommandLine(int argc, char *argv[]) {
    if (Me==RootProc) {
        fprintf(stdout,"+ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +\n");
        fprintf(stdout,"|                   Mechanical AFM Model                      |\n");
        fprintf(stdout,"|  Based on: P. Hapala et al, Phys. Rev. B, 90:085421 (2014)  |\n");
#if SERIAL
        fprintf(stdout,"|           This C implemenation by Peter Spijker             |\n");
#else
        fprintf(stdout,"|         This MPI-C implemenation by Peter Spijker           |\n");
#endif
        fprintf(stdout,"|          2014-2015 (c) Aalto University, Finland            |\n");
        fprintf(stdout,"+ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +\n");
    }
    if ((argc < 2) || (argc > 2)) {
        error("Specify an input file to be read!");
    }
    else {
        sprintf(InputFileName,"%s",argv[1]);
    }
    return;
}

/* A function to read an input file */
void readInputFile(void) {

    FILE *fp;
    int check;
    char keyword[NAME_LENGTH];
    char value[NAME_LENGTH];
    char line[LINE_LENGTH];
    char tmp_coulomb[NAME_LENGTH], tmp_minterm[NAME_LENGTH];
    char tmp_gzip[NAME_LENGTH], tmp_units[NAME_LENGTH];
    char tmp_flexible[NAME_LENGTH], tmp_rigidgrid[NAME_LENGTH];

    /* Initialize the mandatory options */
    sprintf(Options.xyzfile,"");
    sprintf(Options.paramfile,"");
    sprintf(Options.tipatom,"");
    sprintf(Options.dummyatom,"");
    Options.minterm = -1;

    /* Initialize the other options */
    sprintf(Options.planeatom,"");
    Options.units = U_KCAL;
    sprintf(tmp_units,"%s","kcal/mol");
    Options.coulomb = FALSE;
    Options.dx = 0.1;
    Options.dy = 0.1;
    Options.dz = 0.1;
    Options.zlow = 6.0;
    Options.zhigh = 10.0;
    Options.zplane = NEGVAL;
    Options.etol = 0.01;
    Options.ftol = 0.01;
    Options.cfac = 0.001;
    Options.maxsteps = 5000;
    Options.bufsize = 1000;
    Options.gzip = TRUE;
    Options.flexible = FALSE;
    Options.rigidgrid = FALSE;

    /* Check if the file exists */
    fp = fopen(InputFileName,"r");
    if (fp==NULL) {
        error("The file %s does not exist!",InputFileName);
    }

    /* Scan the file line by line */
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Get the keyword and convert to uppercase */
        sscanf(line, "%s %s", keyword, value);
        strlow(keyword);
        if (strcmp(keyword,"xyzfile")==0) {
            sprintf(Options.xyzfile,"%s",value);
        }
        else if (strcmp(keyword,"paramfile")==0) {
            sprintf(Options.paramfile,"%s",value);
        }
        else if (strcmp(keyword,"tipatom")==0) {
            sprintf(Options.tipatom,"%s",value);
        }
        else if (strcmp(keyword,"dummyatom")==0) {
            sprintf(Options.dummyatom,"%s",value);
        }
        else if (strcmp(keyword,"planeatom")==0) {
            sprintf(Options.planeatom,"%s",value);
        }
        else if (strcmp(keyword,"dx")==0) {
            Options.dx = atof(value);
        }
        else if (strcmp(keyword,"dy")==0) {
            Options.dy = atof(value);
        }
        else if (strcmp(keyword,"dz")==0) {
            Options.dz = atof(value);
        }
        else if (strcmp(keyword,"zlow")==0) {
            Options.zlow = atof(value);
        }
        else if (strcmp(keyword,"zhigh")==0) {
            Options.zhigh = atof(value);
        }
        else if (strcmp(keyword,"zplane")==0) {
            Options.zplane = atof(value);
        }
        else if (strcmp(keyword,"etol")==0) {
            Options.etol = atof(value);
        }
        else if (strcmp(keyword,"ftol")==0) {
            Options.ftol = atof(value);
        }
        else if (strcmp(keyword,"cfac")==0) {
            Options.cfac = atof(value);
        }
        else if (strcmp(keyword,"maxsteps")==0) {
            Options.maxsteps = atoi(value);
        }
        else if (strcmp(keyword,"bufsize")==0) {
            Options.bufsize = atoi(value);
        }
        else if (strcmp(keyword,"gzip")==0) {
            if (strcmp(value,"on")==0) {
                Options.gzip = TRUE;
            }
            else if (strcmp(value,"off")==0) {
                Options.gzip = FALSE;
            }
            else {
                error("Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword,"flexible")==0) {
            if (strcmp(value,"on")==0) {
                Options.flexible = TRUE;
            }
            else if (strcmp(value,"off")==0) {
                Options.flexible = FALSE;
            }
            else {
                error("Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword,"rigidgrid")==0) {
            if (strcmp(value,"on")==0) {
                Options.rigidgrid = TRUE;
            }
            else if (strcmp(value,"off")==0) {
                Options.rigidgrid = FALSE;
            }
            else {
                error("Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword,"coulomb")==0) {
            if (strcmp(value,"on")==0) {
                Options.coulomb = TRUE;
            }
            else if (strcmp(value,"off")==0) {
                Options.coulomb = FALSE;
            }
            else {
                error("Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword,"minterm")==0) {
            if (strcmp(value,"e")==0) {
                Options.minterm = MIN_E;
            }
            else if (strcmp(value,"f")==0) {
                Options.minterm = MIN_F;
            }
            else if (strcmp(value,"ef")==0) {
                Options.minterm = MIN_EF;
            }
            else {
                error("Option %s must be either e, f or ef!", keyword);
            }
            sprintf(tmp_minterm,"%s",value);
        }
        else if (strcmp(keyword,"units")==0) {
            if (strcmp(value,"kcal/mol")==0) {
                Options.units = U_KCAL;
            }
            else if (strcmp(value,"kcal")==0) {
                Options.units = U_KCAL;
            }
            else if (strcmp(value,"kJ/mol")==0) {
                Options.units = U_KJ;
            }
            else if (strcmp(value,"kJ")==0) {
                Options.units = U_KJ;
            }
            else if (strcmp(value,"eV")==0) {
                Options.units = U_EV;
            }
            else {
                error("Option %s must be either kcal/mol (default), kJ/mol or eV!", keyword);
            }
            sprintf(tmp_units,"%s",value);
        }
        else {
            error("Unknown option %s!", keyword);
        }
    }

    /* Check if all necessary options are initialized */
    if (strcmp(Options.xyzfile,"")==0) {
        error("Specify at least an xyzfile!");
    }
    if (strcmp(Options.paramfile,"")==0) {
        error("Specify at least a parameter file!");
    }
    if (strcmp(Options.tipatom,"")==0) {
        error("Specify at least a tip atom!");
    }
    if (strcmp(Options.dummyatom,"")==0) {
        error("Specify at least a dummy atom!");
    }
    if (Options.minterm<0) {
        error("Specify at least a minimization termination criterion (e, f, or ef)!");
    }
    if ( (strcmp(Options.planeatom,"")==0) && (Options.zplane <= NEGVAL) ) {
        error("Specify at least a plane or a plane atom!");
    }
    else if ( (strcmp(Options.planeatom,"")!=0) && (Options.zplane > NEGVAL) ) {
        error("Specify only a plane or a plane atom!");
    }

    /* Close file */
    fclose(fp);

    /* Set some useful thingies */
    if (Options.coulomb) {
        sprintf(tmp_coulomb,"%s","on");
    }
    else {
        sprintf(tmp_coulomb,"%s","off");
    }
    if (Options.gzip) {
        sprintf(tmp_gzip,"%s","on");
    }
    else {
        sprintf(tmp_gzip,"%s","off");
    }
    if (Options.flexible) {
        sprintf(tmp_flexible,"%s","on");
    }
    else {
        sprintf(tmp_flexible,"%s","off");
    }
    if (Options.rigidgrid) {
        sprintf(tmp_rigidgrid,"%s","on");
    }
    else {
        sprintf(tmp_rigidgrid,"%s","off");
    }

    /* Do some sanity checking */
    if ((Options.rigidgrid) && (Options.flexible)) {
        error("Cannot use a flexible molecule with a static force grid!");
    }

    /* Set function pointers */
    if (Options.rigidgrid) {
        interactTipSurface = interactTipSurfaceFromGrid;
    }
    else {
        interactTipSurface = interactTipSurfaceDirectly;
    }

    /* Talk to me */
    debugline(RootProc,"");
    debugline(RootProc,"Input settings for %s:", InputFileName);
    debugline(RootProc,"");
    debugline(RootProc,"xyzfile:         %-s", Options.xyzfile);
    debugline(RootProc,"paramfile:   %-s", Options.paramfile);
    debugline(RootProc,"tipatom:         %-s", Options.tipatom);
    debugline(RootProc,"dummyatom:   %-s", Options.dummyatom);
    if (strcmp(Options.planeatom,"")!=0) {
        debugline(RootProc,"planeatom:       %-s", Options.planeatom);
    }
    else if (Options.zplane > NEGVAL) {
        debugline(RootProc,"zplane:          %-8.4f", Options.zplane);
    }
    debugline(RootProc,"");
    debugline(RootProc,"units:           %-s", tmp_units);
    debugline(RootProc,"");
    debugline(RootProc,"minterm:         %-s", tmp_minterm);
    debugline(RootProc,"etol:                %-8.4f", Options.etol);
    debugline(RootProc,"ftol:                %-8.4f", Options.ftol);
    debugline(RootProc,"cfac:                %-8.4f", Options.cfac);
    debugline(RootProc,"maxsteps:        %-8d", Options.maxsteps);
    debugline(RootProc,"");
    debugline(RootProc,"zhigh:           %-8.4f", Options.zhigh);
    debugline(RootProc,"zlow:                %-8.4f", Options.zlow);
    debugline(RootProc,"dx:                  %-8.4f", Options.dx);
    debugline(RootProc,"dy:                  %-8.4f", Options.dy);
    debugline(RootProc,"dz:                  %-8.4f", Options.dz);
    debugline(RootProc,"");
    debugline(RootProc,"coulomb:         %-s",tmp_coulomb);
    debugline(RootProc,"");
    debugline(RootProc,"flexible:        %-s",tmp_flexible);
    debugline(RootProc,"rigidgrid:   %-s",tmp_rigidgrid);
    debugline(RootProc,"");
    debugline(RootProc,"bufsize:         %-8d", Options.bufsize);
    debugline(RootProc,"gzip:                %-s", tmp_gzip);
    debugline(RootProc,"");

    /* Return home */
    return;
}

/* Read the XYZ file */
void readXYZFile(void) {

    FILE *fp;
    int i, n, nplaneatoms, firstline, realxyz, ncols;
    char line[LINE_LENGTH], value[NAME_LENGTH], dump[LINE_LENGTH], *pch;
    double fdump, avgz;

    /* Read the file once, to determine the number of atoms */
    fp = fopen(Options.xyzfile,"r");
    if (fp==NULL) {
        error("No such file: %s!", Options.xyzfile);
    }
    Natoms = 0;
    firstline = TRUE;
    realxyz = FALSE;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* If the first line contains an integer, it is most likely a proper XYZ file */
        if (firstline) {
            sscanf(line, "%s", value);
            if (isint(value)) {
                Natoms = atoi(value);
                realxyz = TRUE;
                break;
            }
            firstline = FALSE;
        }
        /* Otherwise continue reading the lines and count how many atoms we have */
        if (!firstline) {
            /* Count useful lines */
            Natoms++;
        }
    }
    rewind(fp);

    //if (Me == RootProc) fprintf(stdout,"\n*** natoms = %d\n\n",Natoms);

    /* Initialize the global surface atoms vector and lists */
    Surf_pos = (VECTOR *)malloc(Natoms*sizeof(VECTOR));
    Surf_q = (double *)malloc(Natoms*sizeof(double));
    Surf_mass = (double *)malloc(Natoms*sizeof(double));
    Surf_type = (char **)malloc(Natoms*sizeof(char*));
    Surf_fix = (int *)malloc(Natoms*sizeof(int));
    Surf_vel = (VECTOR *)malloc(Natoms*sizeof(VECTOR));
    Surf_force = (VECTOR *)malloc(Natoms*sizeof(VECTOR));
    for (i=0; i<Natoms; ++i) {
        Surf_type[i] = (char *)malloc(ATOM_LENGTH*sizeof(char));
        Surf_q[i] = 0.0;
        Surf_fix[i] = 0;
        Surf_mass[i] = 0.0;
    }

    /* Read each line and store the data */
    i = n = ncols = 0;
    firstline = TRUE;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* If it is a real XYZ file, skip the first two lines */
        if ( (realxyz) && (n<2) ) {
            n++;
            continue;
        }
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Based on the first line with actual atom information, determine how many columns there are */
        if (firstline) {
            strcpy(dump,line);
            pch = strtok(dump," \t\n\r\f");
            while (pch!=NULL) {
                pch = strtok(NULL," \t\n\r\f");
                ncols++;
            }
            firstline = FALSE;
        }
        if (ncols == 4) {
            sscanf(line, "%s %lf %lf %lf", Surf_type[i], &(Surf_pos[i].x), &(Surf_pos[i].y), &(Surf_pos[i].z));
        }
        if (ncols == 5) {
            sscanf(line, "%s %lf %lf %lf %lf", Surf_type[i], &(Surf_pos[i].x), &(Surf_pos[i].y), &(Surf_pos[i].z), &(Surf_q[i]));
        }
        if (ncols == 6) {
            sscanf(line, "%s %lf %lf %lf %lf %d", Surf_type[i], &(Surf_pos[i].x), &(Surf_pos[i].y), &(Surf_pos[i].z), &(Surf_q[i]), &(Surf_fix[i]));
        }
        i++;
    }

    for (i=0; i<Natoms; ++i) {
        if (Surf_fix[i]>0) {
            Nfixed++;
        }
    }

    //if ( Me == RootProc ) {
    //  fprintf(stdout,"\n[%d] *** ncols = %d\n\n",Me,ncols);
    //  for (i=0; i<Natoms; ++i) {
    //  fprintf(stdout,"[%d] >>> %s %lf %lf %lf %lf\n", Me, Surf_type[i], Surf_pos[i].x, Surf_pos[i].y, Surf_pos[i].z, Surf_q[i]);
    //  }
    //}

    /* Put the plane atoms at 0 (in z) or as specified separately*/
    if (strcmp(Options.planeatom,"")!=0) {
        avgz = 0.0;
        nplaneatoms = 0;
        for (i=0; i<Natoms; ++i) {
            if (strcmp(Surf_type[i],Options.planeatom)==0) {
                nplaneatoms++;
                avgz += Surf_pos[i].z;
            }
        }
        avgz /= nplaneatoms;
        for (i=0; i<Natoms; ++i) {
            Surf_pos[i].z -= avgz;
        }
    }
    else if (Options.zplane > NEGVAL) {
        avgz = 0.0;
        for (i=0; i<Natoms; ++i) {
            avgz += Surf_pos[i].z;
        }
        avgz /= Natoms;
        for (i=0; i<Natoms; ++i) {
            Surf_pos[i].z -= avgz;
        }
        for (i=0; i<Natoms; ++i) {
            Surf_pos[i].z += Options.zplane;
        }
    }

    /* Return home */
    return;
}

/* Mixing rule functions */
double mixsig(double sig1, double sig2) {
    return (sig1+sig2)/2;
}
double mixeps(double eps1, double eps2) {
    return sqrt(eps1*eps2);
}

/* Retrieve the index of the specific atom type */
int type2num(char *atom) {
    int i;
    for (i=0; i<Ntypes; ++i) {
        if (strcmp(SurfType2Num[i],atom)==0) {
            break;
        }
    }
    return i;
}

/* Read the parameter file */
void readParameterFile(void) {

    FILE *fp;
    char atom[ATOM_LENGTH], keyword[NAME_LENGTH], dump[NAME_LENGTH], line[LINE_LENGTH];
    char atom1[ATOM_LENGTH], atom2[ATOM_LENGTH], style[NAME_LENGTH], *pch;
    double eps, sig, eps_cross, sig_cross, mass, es6, es12;
    double sig6, De, a, re;
    double eps_tip, sig_tip, q_tip, qbase, epermv;
    int i, j, k, check, hcheck, nplaneatoms, natoms, ncols;
    double avgx, avgy, dx, dy, chargecheck, qdump;

    /* Initialize the universe */
    Box.x = Box.y = Box.z = -1.0;

    /* Open the parameter file */
    fp = fopen(Options.paramfile,"r");
    if (fp==NULL) {
        error("No parameter file %s found!",Options.paramfile);
    }

    /* Scan the parameter file for the universe size and for the tip atom definitions */
    check = FALSE;
    Ntypes = 0;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Read line to determine keyword */
        sscanf(line,"%s",keyword);
        /* Process the separate keywords */
        if (strcmp(keyword,"box")==0) {
            if ( (Box.x<0) && (Box.y<0) && (Box.z<0)) {
                sscanf(line,"%s %lf %lf %lf",dump,&(Box.x),&(Box.y),&(Box.z));
            }
            else {
                error("Keyword box cannot be defined more than once in parameter file!");
            }
        }
        if (strcmp(keyword,"atom")==0) {
            sscanf(line,"%s %s",dump,atom);
            if (strcmp(atom,Options.tipatom)==0) {
                if (check == TRUE) {
                    error("Parameters for tip atom can only be specified once!");
                }
                sscanf(line,"%s %s %lf %lf %s %lf",dump,dump,&(eps_tip),&(sig_tip),dump,&(q_tip));
                check = TRUE;
            }
            Ntypes++;
        }
    }
    if (!check) {
        error("Parameters for tip atom not defined in parameter file!");
    }
    rewind(fp);

    /* Now we know the size of the universe, put the molecule in the center of it */
    if (strcmp(Options.planeatom,"")!=0) {
        avgx = avgy = 0.0;
        nplaneatoms = 0;
        for (i=0; i<Natoms; ++i) {
            if (strcmp(Surf_type[i],Options.planeatom)==0) {
                nplaneatoms++;
                avgx += Surf_pos[i].x;
                avgy += Surf_pos[i].y;
            }
        }
        avgx /= nplaneatoms;
        avgy /= nplaneatoms;
        dx = (Box.x/2) - avgx;
        dy = (Box.y/2) - avgy;
        for (i=0; i<Natoms; ++i) {
            Surf_pos[i].x += dx;
            Surf_pos[i].y += dy;
        }
    }

    /* Set up the interaction list */
    TipSurfParams = (InteractionList *)malloc(Natoms*sizeof(InteractionList));

    /* Create a list with all the possible surface atom particle types */
    Ntypes -= 2;    /* Subtract the tip and dummy atoms */
    if (Ntypes<1) {
        error("Either the tip and/or dummy atom is not specified, or there is no molecule defined. Fix it!");
    }
    SurfSurfParams = (InteractionList **)malloc(Ntypes*sizeof(InteractionList *));
    SurfType2Num = (char **)malloc(Ntypes*sizeof(char *));
    for (i=0; i<Ntypes; ++i) {
        SurfSurfParams[i] = (InteractionList *)malloc(Ntypes*sizeof(InteractionList));
        SurfType2Num[i] = (char *)malloc(ATOM_LENGTH*sizeof(char));
    }
    for (j=0, k=0; j<Natoms; ++j) {
        check = FALSE;
        for (i=0; i<Ntypes; ++i) {
            if (strcmp(Surf_type[j],SurfType2Num[i])==0) {
                check = TRUE;
            }
        }
        if (!check) {
            strcpy(SurfType2Num[k],Surf_type[j]);
            k++;
        }
    }
    if (k!=Ntypes) {
        error("Lost an atom type somewhere along the way");
    }

    /* The constant part of the Coulomb equation (depends on chosen unit system) */
    epermv = 1.0;
    if (Options.units == U_KCAL) {
        qbase = 332.06371;
    }
    else if (Options.units == U_KJ) {
        qbase = 1389.354563;
    }
    else if (Options.units == U_EV) {
        qbase = 14.39964901;
    }
    qbase /= epermv;

    /* Quickly check whether charges were read from the XYZ file                                                     */
    /* NOTE: if only zero charges are specified in the XYZ file, this check sort of fails, */
    /*           as the RMSQE will be zero in that case. If you want zero charge, make sure      */
    /*           the charges in the parameter file are also set to zero!                                             */
    chargecheck = 0.0;
    for (i=0; i<Natoms; ++i) {
        chargecheck += (Surf_q[i]*Surf_q[i]);
    }
    if (fabs(chargecheck)<TOLERANCE) {
        warning("The RMSQ-error for the charges read from the XYZ-file is zero. Charges will be read from the parameter file. If you want zero charge, set the charge in the parameter file to zero.");
    }

    /* Read the parameter file again, but this time, create the interaction list
         for the surface atoms, for the dummy atom, and also for the harmonic spring */
    check = hcheck = FALSE;
    natoms = 0;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Read line to determine keyword */
        sscanf(line,"%s",keyword);
        if (strcmp(keyword,"atom")==0) {
            sscanf(line,"%s %s %lf %lf %lf %lf",dump,atom,&(eps),&(sig),&(mass),&(qdump));
            /* Loop all atoms in the surface and check if they match this parameter set */
            for (i=0; i<Natoms; ++i) {
                if (strcmp(Surf_type[i],atom)==0) {
                    natoms++;
                    eps_cross = mixeps(eps,eps_tip);
                    sig_cross = mixsig(sig,sig_tip);                            /* To power 1 */
                    sig_cross = (sig_cross*sig_cross*sig_cross);    /* To power 3 */
                    sig_cross *= sig_cross;                                             /* To power 6 */
                    TipSurfParams[i].es12 = 4 * eps_cross * sig_cross * sig_cross;
                    TipSurfParams[i].es6    = 4 * eps_cross * sig_cross;
                    if (fabs(chargecheck)<TOLERANCE) {
                        Surf_q[i] = qdump;
                    }
                    TipSurfParams[i].qq     = qbase * q_tip * Surf_q[i];
                    TipSurfParams[i].morse = FALSE;
                    Surf_mass[i] = mass;
                    /* While we read these parameters and assign the tip-molecule interactions,
                        we can also build the molecule-molecule interactions (at least the diagonal) */
                    k = type2num(atom);
                    SurfSurfParams[k][k].eps = eps;
                    SurfSurfParams[k][k].sig = sig;
                }
            }
            /* We found a dummy atom in the parameter list */
            if (strcmp(Options.dummyatom,atom)==0) {
                if (check == TRUE) {
                    error("Parameters for dummy atom can only be specified once!");
                }
                check = TRUE;
                eps_cross = mixeps(eps,eps_tip);
                sig_cross = mixsig(sig,sig_tip);                            /* To power 1 */
                sig_cross = (sig_cross*sig_cross*sig_cross);    /* To power 3 */
                sig_cross *= sig_cross;                                             /* To power 6 */
                DummyParams.es12 = 4 * eps_cross * sig_cross * sig_cross;
                DummyParams.es6  = 4 * eps_cross * sig_cross;
                DummyParams.qq   = 0.0; /* Ignore Coulomb interaction between tip and dummy */
                DummyParams.rmin = mixsig(sig,sig_tip) * SIXTHRT2; /* Needed for tip positioning */
                DummyParams.morse = FALSE;
            }
        }
        if (strcmp(keyword,"harm")==0) {
            if (hcheck == TRUE) {
                error("Parameters for harmonic spring can only be specified once!");
            }
            sscanf(line,"%s %s %lf %lf",dump,atom,&(Harmonic.k),&(Harmonic.r0));
            Harmonic.morse = FALSE;
            if (strcmp(atom,Options.tipatom)!=0) {
                error("Harmonic spring should be defined on tip atom!");
            }
            hcheck = TRUE;
        }
    }
    if (natoms != Natoms) {
        error("Not all atoms have been assigned parameters! (%d/%d)",natoms,Natoms);
    }
    if (check == FALSE) {
        error("Parameters for dummy atom not defined in parameter file!");
    }
    if (hcheck == FALSE) {
        error("No harmonic spring parameters found in parameter file!");
    }

    /* Build the entire interaction matrix for surface-surface interactions */
    for (i=0; i<Ntypes; ++i) {
        for (j=0; j<Ntypes; ++j) {
            eps_cross = mixeps(SurfSurfParams[i][i].eps,SurfSurfParams[j][j].eps);
            sig_cross = mixsig(SurfSurfParams[i][i].sig,SurfSurfParams[j][j].sig);
            // sig_cross = mixsig(sig,sig_tip);                            /* To power 1 */
            sig_cross = (sig_cross*sig_cross*sig_cross);    /* To power 3 */
            sig_cross *= sig_cross;                                             /* To power 6 */
            SurfSurfParams[i][j].es12 = 4 * eps_cross * sig_cross * sig_cross;
            SurfSurfParams[i][j].es6    = 4 * eps_cross * sig_cross;
            SurfSurfParams[j][i].es12 = SurfSurfParams[i][j].es12;
            SurfSurfParams[j][i].es6    = SurfSurfParams[i][j].es6;
            SurfSurfParams[i][j].morse = FALSE;
            SurfSurfParams[j][i].morse = SurfSurfParams[i][j].morse;
        }
    }

    /* Finally check for pair style overwrites */
    rewind(fp);
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Read line to determine keyword */
        sscanf(line,"%s",keyword);
        if (strcmp(keyword,"pair_ovwrt")==0) {
            sscanf(line,"%s %s %s %s",dump,atom1,atom2,style);
            /* Check number of columns */
            ncols = 0;
            strcpy(dump,line);
            pch = strtok(dump," \t\n\r\f");
            while (pch!=NULL) {
                pch = strtok(NULL," \t\n\r\f");
                ncols++;
            }
            /* Lennard-Jones potential */
            if (strcmp(style,"lj")==0) {
                /* Check if number of columns is correct for LJ */
                if (ncols>(4+2)) {
                    error("Only two parameters (eps,sig) allowed for LJ");
                }
                /* Read overwrite parameters */
                sscanf(line,"%s %s %s %s %lf %lf",dump,dump,dump,dump,&(eps),&(sig));
                sig6 = sig*sig*sig;
                sig6 *= sig6;
                es12 = 4 * eps * sig6 * sig6;
                es6  = 4 * eps * sig6;
            }
            /* Morse potential */
            else if (strcmp(style,"morse")==0) {
                /* Check if number of columns is correct for LJ */
                if (ncols>(4+3)) {
                    error("Only three parameters (De,a,re) allowed for Morse");
                }
                /* Read overwrite parameters */
                sscanf(line,"%s %s %s %s %lf %lf %lf",dump,dump,dump,dump,&(De),&(a),&(re));
            }
            /* Try and catch */
            else {
                error("Unknown pair style overwrite '%s'",style);
            }
            /* Store the changes */
            if ((strcmp(atom1,Options.tipatom)==0) || (strcmp(atom2,Options.tipatom)==0)) {
                if (strcmp(atom1,Options.tipatom)==0) {
                    strcpy(atom,atom2);
                }
                if (strcmp(atom2,Options.tipatom)==0) {
                    strcpy(atom,atom1);
                }
                /* Lennard-Jones potential */
                if (strcmp(style,"lj")==0) {
                    for (i=0; i<Natoms; ++i) {
                        if (strcmp(Surf_type[i],atom)==0) {
                            TipSurfParams[i].es12 = es12;
                            TipSurfParams[i].es6    = es6;
                            TipSurfParams[i].morse = FALSE;
                        }
                    }
                }
                /* Morse potential */
                if (strcmp(style,"morse")==0) {
                    for (i=0; i<Natoms; ++i) {
                        if (strcmp(Surf_type[i],atom)==0) {
                            TipSurfParams[i].De = De;
                            TipSurfParams[i].a  = a;
                            TipSurfParams[i].re = re;
                            TipSurfParams[i].morse = TRUE;
                        }
                    }
                }
            }
            else {
                i = type2num(atom1);
                j = type2num(atom2);
                if (strcmp(style,"lj")==0) {
                    SurfSurfParams[i][j].es12 = es12;
                    SurfSurfParams[i][j].es6 = es6;
                    SurfSurfParams[j][i].es12 = SurfSurfParams[i][j].es12;
                    SurfSurfParams[j][i].es6 = SurfSurfParams[i][j].es6;
                    SurfSurfParams[i][j].morse = FALSE;
                    SurfSurfParams[j][i].morse = SurfSurfParams[i][j].morse;
                }
                if (strcmp(style,"morse")==0) {
                    SurfSurfParams[i][j].De = De;
                    SurfSurfParams[i][j].a  = a;
                    SurfSurfParams[i][j].re = re;
                    SurfSurfParams[j][i].De = SurfSurfParams[i][j].De;
                    SurfSurfParams[j][i].a  = SurfSurfParams[i][j].a;
                    SurfSurfParams[j][i].re = SurfSurfParams[i][j].re;
                    SurfSurfParams[i][j].morse = TRUE;
                    SurfSurfParams[j][i].morse = SurfSurfParams[i][j].morse;
                }
            }
        }
    }

    /* Debug printing */
    //for (i=0; i<Ntypes; ++i) {
    //  for (j=0; j<Ntypes; ++j) {
    //      fprintf(stdout,"%d %d - %8.4f %8.4f - %8.4f %8.4f %8.4f\n",i,j,SurfSurfParams[i][j].es12,SurfSurfParams[i][j].es6,SurfSurfParams[i][j].De,SurfSurfParams[i][j].a,SurfSurfParams[i][j].re);
    //  }
    //}

    /* Close file */
    fclose(fp);

    /* Return home */
    return;
}


/***************************
 ** INTERACTION FUNCTIONS **
 ***************************/

/* Interaction between tip and surface */
void interactTipSurfaceDirectly(void) {

    int i;
    double dx, dy, dz, rsqt, dr, dexp, frcfac;
    double e, fx, fy, fz;
    double epair, fpair, terma, termb, termc, sr6, sr12;

    /* Zero the forces and energy */
    TipSurf_energy = 0.0;
    TipSurf_force = NULL_vector;

    /* Zero the forces for the molecule if it's supposed to be flexible */
    if (Options.flexible) {
        for (i=0; i<Natoms; ++i) {
            Surf_force[i].x = 0.0;
            Surf_force[i].y = 0.0;
            Surf_force[i].z = 0.0;
        }
    }

    /* Loop all surface particles */
    for (i=0; i<Natoms; ++i) {
        /* Compute distance (components) */
        dx = Tip_pos.x - Surf_pos[i].x;
        dy = Tip_pos.y - Surf_pos[i].y;
        dz = Tip_pos.z - Surf_pos[i].z;
        rsqt = dx*dx + dy*dy + dz*dz;
        /* Compute the Van der Waals contribution factors */
        if (TipSurfParams[i].morse) {
            /* The LJ is replaced by a Morse potential */
            dr = sqrt(rsqt) - TipSurfParams[i].re;
            dexp = exp(-TipSurfParams[i].a * dr);
            frcfac = 2.0 * TipSurfParams[i].De * TipSurfParams[i].a;
            epair = TipSurfParams[i].De * ( dexp*dexp - 2.0*dexp + 1 );
            fpair = frcfac * ( dexp*dexp - dexp )/sqrt(rsqt);
        }
        else {
            /* The Lennard-Jones 12-6 interaction coefficients */
            sr6 = rsqt*rsqt*rsqt;
            sr12 = sr6*sr6;
            terma = TipSurfParams[i].es12/sr12;
            termb = TipSurfParams[i].es6/sr6;
            epair = terma-termb;
            fpair = (12*terma-6*termb)/rsqt;
        }
        /* The Coulomb interaction coefficients (only if coulomb is on) */
        if (Options.coulomb) {
            termc = TipSurfParams[i].qq/sqrt(rsqt);
        }
        else {
            termc = 0.0;
        }
        /* The interaction energy */
        TipSurf_energy += epair + termc;
        /* The interaction force */
        TipSurf_force.x += (fpair+(termc/rsqt))*dx;
        TipSurf_force.y += (fpair+(termc/rsqt))*dy;
        TipSurf_force.z += (fpair+(termc/rsqt))*dz;
        /* If the molecule is flexible, also store its forces */
        if (Options.flexible) {
            Surf_force[i].x -= (fpair+(termc/rsqt))*dx;
            Surf_force[i].y -= (fpair+(termc/rsqt))*dy;
            Surf_force[i].z -= (fpair+(termc/rsqt))*dz;
        }
    }

    /* Go home */
    return;
}

/* Interaction between tip and dummy atom */
void interactTipDummy(void) {

    double dx, dy, dz, rsqt, dr, dexp, frcfac;
    double epair, fpair, terma, termb, termc, sr6, sr12;

    /* Compute distance (components) */
    dx = Tip_pos.x - Dummy_pos.x;
    dy = Tip_pos.y - Dummy_pos.y;
    dz = Tip_pos.z - Dummy_pos.z;
    rsqt = dx*dx + dy*dy + dz*dz;
    /* The Lennard-Jones interaction coefficients */
    sr6 = rsqt*rsqt*rsqt;
    sr12 = sr6*sr6;
    terma = DummyParams.es12/sr12;
    termb = DummyParams.es6/sr6;
    epair = terma-termb;
    fpair = (12*terma-6*termb)/rsqt;
    /* Now check if the LJ is replaced by a Morse potential */
    if (DummyParams.morse) {
        dr = sqrt(rsqt) - DummyParams.re;
        dexp = exp(-DummyParams.a * dr);
        frcfac = 2.0 * DummyParams.De * DummyParams.a;
        epair = DummyParams.De * ( dexp*dexp - 2.0*dexp + 1 );
        fpair = frcfac * ( dexp*dexp - dexp )/sqrt(rsqt);
    }
    /* The interaction energy */
    TipDummy_energy = epair;
    /* The interaction force */
    TipDummy_force.x = fpair*dx;
    TipDummy_force.y = fpair*dy;
    TipDummy_force.z = fpair*dz;

    /* Go home */
    return;
}

/* Interaction between tip and harmonic constraint */
void interactTipHarmonic(void) {

    double dx, dy, r, dr, rk, fharm;

    /* Compute distance (components, but not in z!) */
    dx = Tip_pos.x - Dummy_pos.x;
    dy = Tip_pos.y - Dummy_pos.y;
    r = sqrt(dx*dx + dy*dy);
    /* Compute the harmonic coefficients */
    dr = r - Harmonic.r0;
    rk = Harmonic.k*dr;
    /* The interaction energy */
    TipHarmonic_energy = rk*dr;
    /* The interaction force */
    if (r>TOLERANCE) {
        fharm = (-2*rk/r);
    }
    else {
        fharm = 0.0;
    }
    TipHarmonic_force.x = fharm*dx;
    TipHarmonic_force.y = fharm*dy;
    TipHarmonic_force.z = 0.0;

    /* Go home */
    return;
}

/* Build a 3D force grid (for efficient look ups in case there is no flexibility) */
void build3DForceGrid(void) {

    double mingridspacing, smallestgridspacing;
    double *localforcegrid;
    int ix, iy, iz;
    int i, n, onproc;
    int nx, ny, nz;
    double x, y, z;
    double zoffset;
    FILE *fp;

    /* Check the grid spacing, and adjust if necessary */
    mingridspacing = 9e99;
    smallestgridspacing = 0.01; /* A grid smaller than this (in Angstrom) is silly */
    if (Options.dx < mingridspacing) {
        mingridspacing = Options.dx;
    }
    if (Options.dy < mingridspacing) {
        mingridspacing = Options.dy;
    }
    if (Options.dz < mingridspacing) {
        mingridspacing = Options.dz;
    }
    if (mingridspacing < smallestgridspacing) {
        mingridspacing = smallestgridspacing;
    }
    GridSpacing = mingridspacing;

    /* How many points in the grid in each direction */
    Ngrid.x = (int)((Npoints.x*Options.dx)/GridSpacing);
    Ngrid.y = (int)((Npoints.y*Options.dy)/GridSpacing);
    Ngrid.z = (int)((Npoints.z*Options.dz)/GridSpacing);
    Ngridpoints = (Ngrid.x+1)*(Ngrid.y+1)*(Ngrid.z+1);

    /* Allocate memory for force grids (locally and globally) */
    ForceGridRigid = (double *)malloc(3*Ngridpoints*sizeof(double));
    localforcegrid = (double *)malloc(3*Ngridpoints*sizeof(double));
    for (i=0; i<(3*Ngridpoints); ++i) {
        ForceGridRigid[i] = 0.0;
        localforcegrid[i] = 0.0;
    }

    /* Loop all possible grid points and locally compute the VdW force */
    debugline(RootProc,"Computing 3D force grid (%d grid points)",Ngridpoints);
    zoffset = Options.zhigh - DummyParams.rmin - Ngrid.z*GridSpacing;
    /* First over x */
    for (ix=0; ix<=Ngrid.x; ++ix) {
        /* Where is x */
        Tip_pos.x = ix*GridSpacing;
        /* Then over y */
        for (iy=0; iy<=Ngrid.y; ++iy) {
            /* Where is y */
            Tip_pos.y = iy*GridSpacing;
            /* Which processor are we on? Do we skip this point? */
            n = ix*(Ngrid.y+1) + iy;
            onproc = n % NProcessors;
            if (onproc != Me) {
                continue;
            }
            /* Finally loop over z */
            for (iz=0; iz<=Ngrid.z; ++iz) {
                /* Where is z */
                Tip_pos.z = iz*GridSpacing + zoffset;
                /* Compute the VdW/Coulomb interaction */
                /* We cannot use the function pointer here, this has to be explicit! */
                interactTipSurfaceDirectly();
                /* Which point in the force grid are we storing */
                nx = ((ix*(Ngrid.y+1)) + iy)*(Ngrid.z+1) + iz;
                ny = nx + Ngridpoints;
                nz = ny + Ngridpoints;
                /* Store the force */
                localforcegrid[nx] = TipSurf_force.x;
                localforcegrid[ny] = TipSurf_force.y;
                localforcegrid[nz] = TipSurf_force.z;
            } /* End loop in z */
        } /* End loop in y */
    } /* End loop in x */

    /* First wait and then communicate all the data to all processors */
#if !SERIAL
    MPI_Barrier(Universe);
    MPI_Allreduce(localforcegrid,ForceGridRigid,3*Ngridpoints,MPI_DOUBLE,MPI_SUM,Universe);
#else
    for (i=0; i<(3*Ngridpoints); ++i) {
        ForceGridRigid[i] = localforcegrid[i];
    }
#endif

    /* We probably do not need to write the force grid to file. It will be huge and
         because the grid is computer in parallel, it is relatively fast to recompute it. */

    /* The main processor will write the grid to file */
    //if (Me == RootProc) {
    //  fp = fopen("frcgrid.dat","w");
    //  for (ix=0; ix<=Ngrid.x; ++ix) {
    //      for (iy=0; iy<=Ngrid.y; ++iy) {
    //  for (iz=0; iz<=Ngrid.z; ++iz) {
    //          /* Which points do we need to retrieve? */
    //          nx = ((ix*(Ngrid.y+1)) + iy)*(Ngrid.z+1) + iz;
    //      ny = nx + Ngridpoints;
    //      nz = ny + Ngridpoints;
    //      fprintf(fp,"%d %d %d %10.6f %10.6f %10.6f\n",ix,iy,iz,ForceGridRigid[nx],ForceGridRigid[ny],ForceGridRigid[nz]);
    //  }
    //  }
    //  }
    //  fclose(fp);
    //}

    /* Go home */
    return;
}

/* Retrieve grid point ... */
void retrieveGridPoint(IVECTOR *gp,VECTOR *dxyz) {
    double zoffset;

    /* What is the nearest grid point to the current position of the tip */
    zoffset = Options.zhigh - DummyParams.rmin - Ngrid.z*GridSpacing;
    gp->x = (int)(Tip_pos.x/GridSpacing);
    gp->y = (int)(Tip_pos.y/GridSpacing);
    gp->z = (int)((Tip_pos.z-zoffset)/GridSpacing);

    /* Fix possible out of bounds (only happens at the edges) */
    if (Tip_pos.x<0) {
        gp->x = 0;
    }
    if (Tip_pos.y<0) {
        gp->y = 0;
    }
    if (Tip_pos.z<zoffset) {
        gp->z = 0;
    }
    if (Tip_pos.x>Box.x) {
        gp->x = Ngrid.x-1;
    }
    if (Tip_pos.y>Box.y) {
        gp->y = Ngrid.y-1;
    }
    if (Tip_pos.z>(zoffset+Ngrid.z*GridSpacing)) {
        gp->z = Ngrid.z-1;
    }

    /* And how far are we inside the cube (fractional coordinates) */
    dxyz->x = (Tip_pos.x - ((gp->x)*GridSpacing)) / GridSpacing;
    dxyz->y = (Tip_pos.y - ((gp->y)*GridSpacing)) / GridSpacing;
    dxyz->z = ((Tip_pos.z-zoffset) - ((gp->z)*GridSpacing)) / GridSpacing;

    /* Fix possible out of bounds (only happens at the edges) */
    if (Tip_pos.x<0) {
        dxyz->x = 0.0;
    }
    if (Tip_pos.y<0) {
        dxyz->y = 0.0;
    }
    if (Tip_pos.z<zoffset) {
        dxyz->z = 0.0;
    }
    if (Tip_pos.x>Box.x) {
        dxyz->x = 1.0;
    }
    if (Tip_pos.y>Box.y) {
        dxyz->y = 1.0;
    }
    if (Tip_pos.z>(zoffset+Ngrid.z*GridSpacing)) {
        dxyz->z = 1.0;
    }

    /* Go home */
    return;
}

/* ... and corresponding array indices */
void retrieveArrayIndex(IVECTOR gp, IVECTOR *ai) {

    /* Compute the array index based on the given grid point */
    ai->x = ((gp.x*(Ngrid.y+1)) + gp.y)*(Ngrid.z+1) + gp.z;
    ai->y = ai->x + Ngridpoints;
    ai->z = ai->y + Ngridpoints;

    /* Go home */
    return;
}

/* Interaction between tip and surface (interpolation from a grid) */
void interactTipSurfaceFromGrid(void) {

    IVECTOR gp, cgp;
    VECTOR dxyz;
    IVECTOR ai;
    int ix, iy, iz, n;
    double zoffset;
    double *fgr;
    VECTOR C[8];
    double c00, c10, c01, c11, c0, c1, c;

    /* Zero the forces and energy */
    TipSurf_energy = 0.0;
    TipSurf_force = NULL_vector;

    /* What is the nearest grid point to the current position of the tip? */
    retrieveGridPoint(&gp,&dxyz);

    /* Find the surrounding grid points as array indices and retrieve the force values */
    fgr = ForceGridRigid;
    for (ix=0, n=0; ix<=1; ++ix) {
        cgp.x = gp.x + ix;
        for (iy=0; iy<=1; ++iy) {
            cgp.y = gp.y + iy;
            for (iz=0; iz<=1; ++iz) {
                cgp.z = gp.z + iz;
                /* Get the array index */
                retrieveArrayIndex(cgp,&ai);
                /* Retrieve and store the force */
                C[n].x = fgr[ai.x];
                C[n].y = fgr[ai.y];
                C[n].z = fgr[ai.z];
                n++;
            }
        }
    }

    /* This is the order of the coefficients in the C array */
    /*      000 - 001 - 010 - 011 - 100 - 101 - 110 - 111           */

    /* The trilinear interpolation below is based on:                */
    /*  http://en.wikipedia.org/wiki/Trilinear_interpolation */

    /* Construct the force on the tip (in x) */
    c00 = C[0].x*(1-dxyz.x) + C[4].x*(dxyz.x);
    c01 = C[1].x*(1-dxyz.x) + C[5].x*(dxyz.x);
    c10 = C[2].x*(1-dxyz.x) + C[6].x*(dxyz.x);
    c11 = C[3].x*(1-dxyz.x) + C[7].x*(dxyz.x);
    c0 = c00*(1-dxyz.y) + c10*(dxyz.y);
    c1 = c01*(1-dxyz.y) + c11*(dxyz.y);
    TipSurf_force.x = c0*(1-dxyz.z) + c1*(dxyz.z);

    /* Construct the force on the tip (in y) */
    c00 = C[0].y*(1-dxyz.x) + C[4].y*(dxyz.x);
    c01 = C[1].y*(1-dxyz.x) + C[5].y*(dxyz.x);
    c10 = C[2].y*(1-dxyz.x) + C[6].y*(dxyz.x);
    c11 = C[3].y*(1-dxyz.x) + C[7].y*(dxyz.x);
    c0 = c00*(1-dxyz.y) + c10*(dxyz.y);
    c1 = c01*(1-dxyz.y) + c11*(dxyz.y);
    TipSurf_force.y = c0*(1-dxyz.z) + c1*(dxyz.z);

    /* Construct the force on the tip (in z) */
    c00 = C[0].z*(1-dxyz.x) + C[4].z*(dxyz.x);
    c01 = C[1].z*(1-dxyz.x) + C[5].z*(dxyz.x);
    c10 = C[2].z*(1-dxyz.x) + C[6].z*(dxyz.x);
    c11 = C[3].z*(1-dxyz.x) + C[7].z*(dxyz.x);
    c0 = c00*(1-dxyz.y) + c10*(dxyz.y);
    c1 = c01*(1-dxyz.y) + c11*(dxyz.y);
    TipSurf_force.z = c0*(1-dxyz.z) + c1*(dxyz.z);

    /* Go home */
    return;
}

/*********************************
 ** FLEXIBLE MOLECULE FUNCTIONS **
 *********************************/

/* In the case we have a flexible molecule we have to build the topology */
void buildTopology(void) {

    FILE *fp;
    int i, j, k, n, bcheck, acheck, scheck, tcheck, ntopobonds, idmin;
    int a11, a12, a21, a22;
    char keyword[NAME_LENGTH], dump[NAME_LENGTH], line[LINE_LENGTH];
    double kbond, kangle, rzero, d, dx, dy, dz, safedist;
    double theta, dx1, dx2, dy1, dy2, dz1, dz2, d1, d2, costheta;
    double *xrange, *yrange, minx, miny, maxx, maxy, minz;
    char atom1[NAME_LENGTH], atom2[NAME_LENGTH];
    PossibleBonds *posbonds;
    BondInteraction *tmpbonds;
    AngleInteraction *tmpangles;

    /* Open the parameter file */
    fp = fopen(Options.paramfile,"r");

    /* Quickly scan the parameter file to know how many "topobond" keywords there are */
    ntopobonds = 0;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Read line to determine keyword */
        sscanf(line,"%s",keyword);
        /* If it is a topobond, count it */
        if (strcmp(keyword,"topobond")==0) {
            ntopobonds++;
        }
    }
    if (ntopobonds == 0) {
        warning("No topobond keyword found in parameter file, no molecule? Nothing to keep flexible.");
    }
    rewind(fp);

    /* Create an array to keep the possible bonds listed in the parameter file */
    posbonds = (PossibleBonds *)malloc(ntopobonds*sizeof(PossibleBonds));

    /* Read the parameter file again, but this time, parse
         everything we need for modeling a flexible molecule */
    bcheck = acheck = scheck = FALSE;
    ntopobonds = 0;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Read line to determine keyword */
        sscanf(line,"%s",keyword);
        /* The strength of the harmonic bond from the parameter file */
        if (strcmp(keyword,"bond")==0) {
            if (bcheck == TRUE) {
                error("Parameters for harmonic bonds can only be specified once!");
            }
            sscanf(line,"%s %lf",dump,&(kbond));
            bcheck = TRUE;
        }
        /* The strength of the harmonic angle from the parameter file */
        if (strcmp(keyword,"angle")==0) {
            if (acheck == TRUE) {
                error("Parameters for harmonic angles can only be specified once!");
            }
            sscanf(line,"%s %lf",dump,&(kangle));
            acheck = TRUE;
        }
        /* The strength of the harmonic substrate support from the parameter file */
        if (strcmp(keyword,"substrate")==0) {
            if (scheck == TRUE) {
                error("Parameters for harmonic substrate support can only be specified once!");
            }
            sscanf(line,"%s %lf",dump,&(Substrate.k));
            scheck = TRUE;
        }
        /* Collect the bonds, possibly present in the system */
        if (strcmp(keyword,"topobond")==0) {
            tcheck = FALSE;
            sscanf(line,"%s %s %s %lf",dump,atom1,atom2,&(rzero));
            for (i=0; i<ntopobonds; ++i) {
                if (((strcmp(atom1,posbonds[i].a1)==0) && (strcmp(atom2,posbonds[i].a2)==0)) ||
                   ((strcmp(atom1,posbonds[i].a2)==0) && (strcmp(atom2,posbonds[i].a1)==0))) {
                    tcheck = TRUE;
                }
            }
            if (tcheck) {
                error("The topobond for %s and %s is defined at least twice!",atom1,atom2);
            }
            strcpy(posbonds[ntopobonds].a1,atom1);
            strcpy(posbonds[ntopobonds].a2,atom2);
            posbonds[ntopobonds].r0 = rzero;
            ntopobonds++;
        }
    }
    if (bcheck == FALSE) {
        error("No harmonic bond parameters found in parameter file!");
    }
    if (acheck == FALSE) {
        error("No harmonic angle parameters found in parameter file!");
    }
    if (scheck == FALSE) {
        error("No harmonic substrate support parameters found in parameter file!");
    }

    /* Close file */
    fclose(fp);

    /* Some temporary bond storage array */
    tmpbonds = (BondInteraction *)malloc((Natoms*10)*sizeof(BondInteraction));

    /* Figure out whether we need to store this bond */
    safedist = 1.10;
    Nbonds = 0;
    for (i=0; i<Natoms; ++i) {
        for (j=(i+1); j<Natoms; ++j) {
            /* Are the two atoms a possible bond? */
            tcheck = FALSE;
            strcpy(atom1,Surf_type[i]);
            strcpy(atom2,Surf_type[j]);
            for (k=0; k<ntopobonds; ++k) {
                if (((strcmp(atom1,posbonds[k].a1)==0) && (strcmp(atom2,posbonds[k].a2)==0)) ||
                ((strcmp(atom1,posbonds[k].a2)==0) && (strcmp(atom2,posbonds[k].a1)==0))) {
                    tcheck = TRUE;
                    rzero = posbonds[k].r0;
                }
            }
            //fprintf(stdout,"*** %2d %2d - %s %s - %d - %f\n",i,j,atom1,atom2,tcheck,rzero);
            if (tcheck==FALSE) {
                continue;
            }
            /* Compute the distance between the two atoms */
            dx = (Surf_pos[i].x - Surf_pos[j].x);
            dy = (Surf_pos[i].y - Surf_pos[j].y);
            dz = (Surf_pos[i].z - Surf_pos[j].z);
            d = sqrt(dx*dx + dy*dy + dz*dz);
            /* Are these two atoms forming a bond? */
            if (d<(rzero*safedist)) {
                tmpbonds[Nbonds].a1 = i;
                tmpbonds[Nbonds].a2 = j;
                tmpbonds[Nbonds].r0 = d;            /* Yes, we keep the current length as rzero! */
                tmpbonds[Nbonds].k  = kbond;    /* From the parameter file */
                Nbonds++;
            }
        }
    }

    /* Copy all found bonds to a global array */
    Bonds = (BondInteraction *)malloc(Nbonds*sizeof(BondInteraction));
    for (i=0; i<Nbonds; ++i) {
        Bonds[i].a1 = tmpbonds[i].a1;
        Bonds[i].a2 = tmpbonds[i].a2;
        Bonds[i].r0 = tmpbonds[i].r0;
        Bonds[i].k  = tmpbonds[i].k;
        //fprintf(stdout,">>>BOND<<< (%3d/%3d): %3d %3d %f %f\n",i,Nbonds,Bonds[i].a1,Bonds[i].a2,Bonds[i].r0,Bonds[i].k);
    }

    /* Some temporary angle storage array */
    tmpangles = (AngleInteraction *)malloc((Nbonds*10)*sizeof(AngleInteraction));

    /* Find the angles by looping through the bond list */
    Nangles = 0;
    /* Loop all bonds once */
    for (i=0; i<Nbonds; ++i) {
        /* Which atoms form the i-th bond */
        a11 = Bonds[i].a1;
        a12 = Bonds[i].a2;
        /* And loop all bonds twice */
        for (j=(i+1); j<Nbonds; ++j) {
            /* And which atoms form the j-th bond */
            a21 = Bonds[j].a1;
            a22 = Bonds[j].a2;
            /* If both bonds share one and only one atom, an angle is possible */
            /* But only if the other atom is not the same in both bonds */
            if ((a11==a21) && (a12!=a22)) {
                tmpangles[Nangles].a1 = a12;
                tmpangles[Nangles].a2 = a11;
                tmpangles[Nangles].a3 = a22;
            }
            else if ((a11==a22) && (a12!=a21)) {
                tmpangles[Nangles].a1 = a12;
                tmpangles[Nangles].a2 = a11;
                tmpangles[Nangles].a3 = a21;
            }
            else if ((a12==a21) && (a11!=a22)) {
                tmpangles[Nangles].a1 = a11;
                tmpangles[Nangles].a2 = a12;
                tmpangles[Nangles].a3 = a22;
            }
            else if ((a12==a22) && (a11!=a21)) {
                tmpangles[Nangles].a1 = a11;
                tmpangles[Nangles].a2 = a12;
                tmpangles[Nangles].a3 = a21;
            }
            else {
                continue;
            }
            /* We won't get here, unless we found an angle, so store it properly */
            tmpangles[Nangles].bond1 = i;
            tmpangles[Nangles].bond2 = j;
            tmpangles[Nangles].k = kangle;
            /* And compute the base angle */
            dx1 = Surf_pos[tmpangles[Nangles].a1].x - Surf_pos[tmpangles[Nangles].a2].x;
            dx2 = Surf_pos[tmpangles[Nangles].a3].x - Surf_pos[tmpangles[Nangles].a2].x;
            dy1 = Surf_pos[tmpangles[Nangles].a1].y - Surf_pos[tmpangles[Nangles].a2].y;
            dy2 = Surf_pos[tmpangles[Nangles].a3].y - Surf_pos[tmpangles[Nangles].a2].y;
            dz1 = Surf_pos[tmpangles[Nangles].a1].z - Surf_pos[tmpangles[Nangles].a2].z;
            dz2 = Surf_pos[tmpangles[Nangles].a3].z - Surf_pos[tmpangles[Nangles].a2].z;
            d1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
            d2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
            costheta = (dx1*dx2 + dy1*dy2 + dz1*dz2) / (d1*d2);
            if ( costheta >  1.0 ) {
                costheta =  1.0;
            }
            if ( costheta < -1.0 ) {
                costheta = -1.0;
            }
            tmpangles[Nangles].theta0 = acos(costheta);
            /* Increment the angle number */
            Nangles++;
        }
    }

    /* Copy all found angles to a global array */
    Angles = (AngleInteraction *)malloc(Nangles*sizeof(AngleInteraction));
    for (i=0; i<Nangles; ++i) {
        Angles[i].a1 = tmpangles[i].a1;
        Angles[i].a2 = tmpangles[i].a2;
        Angles[i].a3 = tmpangles[i].a3;
        Angles[i].k  = tmpangles[i].k;
        Angles[i].theta0 = tmpangles[i].theta0;
        Angles[i].bond1  = tmpangles[i].bond1;
        Angles[i].bond2  = tmpangles[i].bond2;
        //fprintf(stdout,">>>ANGLE<<< (%3d/%3d): %3d %3d %3d - ",i,Nangles,Angles[i].a1,Angles[i].a2,Angles[i].a3);
        //fprintf(stdout,"%s %s %s - ",Surf_type[Angles[i].a1],Surf_type[Angles[i].a2],Surf_type[Angles[i].a3]);
        //fprintf(stdout,"%f %f (%3d %3d)\n",Angles[i].theta0,Angles[i].k,Angles[i].bond1,Angles[i].bond2);
    }

    /* In case we have a flexible molecule, make a copy of the current surface atom positions */
    if (Options.flexible) {
        Surf_pos_org = (VECTOR *)malloc(Natoms*sizeof(VECTOR));
        for (i=0; i<Natoms; ++i) {
            Surf_pos_org[i].x = Surf_pos[i].x;
            Surf_pos_org[i].y = Surf_pos[i].y;
            Surf_pos_org[i].z = Surf_pos[i].z;
        }
    }

    /* Free some memory */
    free(tmpbonds);
    free(tmpangles);

    /* Return home */
    return;
}

/* Compute the bond contribution within the molecule */
void computeBondsFlexMol(void) {

    int i, a1, a2;
    double f, k, r0, dr, r, dx, dy, dz;

    /* Compute the forces from all internal bonds */
    for (i=0; i<Nbonds; ++i) {
        /* Retrieve the information */
        a1 = Bonds[i].a1;
        a2 = Bonds[i].a2;
        r0 = Bonds[i].r0;
        k  = Bonds[i].k;
        /* Compute the current bond length */
        dx = Surf_pos[a1].x - Surf_pos[a2].x;
        dy = Surf_pos[a1].y - Surf_pos[a2].y;
        dz = Surf_pos[a1].z - Surf_pos[a2].z;
        r = sqrt(dx*dx + dy*dy + dz*dz);
        /* And the harmonic force arising from it */
        dr = r - r0;
        f = -2.0*k*dr/r;
        /* Assign to both atoms */
        Surf_force[a1].x += dx*f;
        Surf_force[a1].y += dy*f;
        Surf_force[a1].z += dz*f;
        Surf_force[a2].x -= dx*f;
        Surf_force[a2].y -= dy*f;
        Surf_force[a2].z -= dz*f;
    }

    /* Go home */
    return;
}

/* Compute the angle contribution within the molecule */
void computeAnglesFlexMol(void) {

    int i, a1, a2, a3;
    double dx1, dx2, dy1, dy2, dz1, dz2, d1, d2;
    double costheta, sintheta, dtheta, theta0;
    double k, tk, f, f11, f12, f22;
    VECTOR f1, f3;

    /* Compute the forces from all internal angles */
    for (i=0; i<Nangles; ++i) {
        /* Retrieve information */
        a1 = Angles[i].a1;
        a2 = Angles[i].a2;
        a3 = Angles[i].a3;
        theta0 = Angles[i].theta0;
        k = Angles[i].k;
        /* Compute the current "angle" */
        dx1 = Surf_pos[a1].x - Surf_pos[a2].x;
        dx2 = Surf_pos[a3].x - Surf_pos[a2].x;
        dy1 = Surf_pos[a1].y - Surf_pos[a2].y;
        dy2 = Surf_pos[a3].y - Surf_pos[a2].y;
        dz1 = Surf_pos[a1].z - Surf_pos[a2].z;
        dz2 = Surf_pos[a3].z - Surf_pos[a2].z;
        d1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        d2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
        costheta = (dx1*dx2 + dy1*dy2 + dz1*dz2) / (d1*d2);
        if ( costheta >  1.0 ) {
            costheta =  1.0;
        }
        if ( costheta < -1.0 ) {
            costheta = -1.0;
        }
        sintheta = sqrt(1.0 - costheta*costheta);
        if ( sintheta < 0.001 ) {
            sintheta = 0.001;
        }
        sintheta = 1.0 / sintheta;
        /* Compute the force */
        dtheta = acos(costheta) - theta0;
        tk = k*dtheta;
        f = -2.0 * tk * sintheta;
        f11 = f*costheta / (d1*d1);
        f12 = -f / (d1*d2);
        f22 = f*costheta / (d2*d2);
        f1.x = f11*dx1 + f12*dx2;
        f1.y = f11*dy1 + f12*dy2;
        f1.z = f11*dz1 + f12*dz2;
        f3.x = f22*dx2 + f12*dx1;
        f3.y = f22*dy2 + f12*dy1;
        f3.z = f22*dz2 + f12*dz1;
        /* Assign to all atoms */
        Surf_force[a1].x += f1.x;
        Surf_force[a1].y += f1.y;
        Surf_force[a1].z += f1.z;
        Surf_force[a2].x -= f1.x + f3.x;
        Surf_force[a2].y -= f1.y + f3.y;
        Surf_force[a2].z -= f1.z + f3.z;
        Surf_force[a3].x += f3.x;
        Surf_force[a3].y += f3.y;
        Surf_force[a3].z += f3.z;
    }

    /* Go home */
    return;

}

/* Compute the contribution of the surface support */
void computeSubstrateSupport(void) {

    int i;
    double dz;

    /* Loop all surface atoms and compute the z-offset with respect to the original coordinates */
    for (i=0; i<Natoms; ++i) {
        /* Deviation from original z-coordinate */
        dz = Surf_pos_org[i].z - Surf_pos[i].z;
        /* If the atom has moved away from the substrate, do nothing */
        if (dz<0) {
            continue;
        }
        /* Assign harmonic force (only in z) */
        Surf_force[i].z += 2.0 * Substrate.k * dz;
    }

    /* BUT WHAT IF THE MOLECULE BOUNCES AWAY ??? */

    /* Return home */
    return;
}

/* Update the positions of the molecules based on the forces, steepest descent minimization */
void updateFlexibleMolecule(void) {

    int i;

    /* Call all force computations for the flexible molecule */
    computeBondsFlexMol();
    computeAnglesFlexMol();
    computeSubstrateSupport();

    /* Loop all surface atoms and update their positions */
    for (i=0; i<Natoms; ++i) {
        Surf_pos[i].x += Options.cfac * Surf_force[i].x;
        Surf_pos[i].y += Options.cfac * Surf_force[i].y;
        Surf_pos[i].z += Options.cfac * Surf_force[i].z;
    }

    /* Keep some atoms fixed even if the molecule is flexible (substrate support) */
    for (i=0; i<Natoms; ++i) {
        if (Surf_fix[i]) {
            Surf_pos[i].x = Surf_pos_org[i].x;
            Surf_pos[i].y = Surf_pos_org[i].y;
            Surf_pos[i].z = Surf_pos_org[i].z;
        }
    }

    /* Return home */
    return;
}


/***************************
 ** FILE OUTPUT FUNCTIONS **
 ***************************/

void dumpToFiles(BUFFER *sendbuf, BUFFER *recvbuf, int bufsize) {

    int i, f, nsr, *curbufsize, *lcbs;
#if !SERIAL
    MPI_Status mpistatus;
#endif

    /* Build an array of the current buffer size for broadcast */
    curbufsize = (int *)malloc(NProcessors*sizeof(int));
    lcbs = (int *)malloc(NProcessors*sizeof(int));
    for (i=0; i<NProcessors; ++i) {
        lcbs[i] = 0;
    }
    lcbs[Me] = bufsize;
#if !SERIAL
    MPI_Allreduce(lcbs,curbufsize,NProcessors,MPI_INT,MPI_SUM,Universe);
#else
    curbufsize[Me] = lcbs[Me];
#endif

#if !SERIAL
    /* Send the data to the root processor */
    if (Me != RootProc) {
        MPI_Send(sendbuf,curbufsize[Me]*sizeof(BUFFER),MPI_CHAR,RootProc,0,Universe);
    }
    /* Receive the date from the daughter processors and write to file */
    else {
#endif
        /* Loop the processors */
        for (i=0; i<NProcessors; ++i) {
            /* On the main processor we have to copy the data only, no send and receive */
            if (i==0) {
                recvbuf = sendbuf;
            }
#if !SERIAL
            /* For all other processors we need to receive the data */
            else {
                MPI_Recv(recvbuf,curbufsize[i]*sizeof(BUFFER),MPI_CHAR,i,0,Universe,&mpistatus);
            }
#endif
            /* Write data to file (only the root processor can do this) */
            /* PLEASE NOTE: DATA IS SENT IN STRIPED FORM, THEY ARE NOT ORDERED! */
            for (nsr=0; nsr<curbufsize[i]; ++nsr) {
                f = recvbuf[nsr].iz;
                /* The file buffer can be a gzip pipe or an ASCII file stream */
                fprintf(FStreams[f],"%d %d %d ",recvbuf[nsr].iz,recvbuf[nsr].ix,recvbuf[nsr].iy);
                fprintf(FStreams[f],"%6.3f %6.3f %6.3f ",recvbuf[nsr].pos.x,recvbuf[nsr].pos.y,recvbuf[nsr].pos.z);
                fprintf(FStreams[f],"%8.4f %8.4f %8.4f ",recvbuf[nsr].f.x,recvbuf[nsr].f.y,recvbuf[nsr].f.z);
                fprintf(FStreams[f],"%6.3f %6.3f %6.3f ",recvbuf[nsr].d.x,recvbuf[nsr].d.y,recvbuf[nsr].d.z);
                fprintf(FStreams[f],"%6.3f %8.4f ",recvbuf[nsr].dd,recvbuf[nsr].angle);
                fprintf(FStreams[f],"%8.4f %d\n",recvbuf[nsr].e,recvbuf[nsr].n);
            }
        }
#if !SERIAL
    }
#endif

    /* Get rid of the buffer size broadcast arrays */
    free(curbufsize);
    free(lcbs);

    /* Go home */
    return;
}

/************************************
 ** FREQUENCY SHIFT APPROXIMATIONS **
 ************************************/

/* void computeDeltaF(double x, double y, VECTOR *ftip) { */

/*   int i; */
/*   double z; */

/*   /\* Wait! *\/ */
/*   MPI_Barrier(Universe); */

/*   /\* Copy data *\/ */
/*   for (i=0; i<=Npoints.z; ++i) { */
/*       z = Options.zhigh - i*Options.dz; */
/*       fprintf(stdout,"@ %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",x,y,z,ftip[i].x,ftip[i].y,ftip[i].z); */
/*   } */

/*   /\* Go home *\/ */
/*   return; */
/* } */

/*************************************************
 ** EVERYTHING TO DO WITH THE SIMULATION ITSELF **
 *************************************************/

/* Set some global variables we need and open the file streams */
void openUniverse(void) {

    int i, n;
    double z;
    char outfile[NAME_LENGTH];

    /* How many points to compute the tip relaxation on */
    Npoints.x = (int)(Box.x/Options.dx);
    Npoints.y = (int)(Box.y/Options.dy);
    Npoints.z = (int)((Options.zhigh-Options.zlow)/Options.dz);
    n = (Npoints.x + 1) * (Npoints.y + 1) * (Npoints.z + 1);
    debugline(RootProc,"3D data grid is: %d x %d x %d (%d in total)",1+Npoints.x,1+Npoints.y,1+Npoints.z,n);

#if !SERIAL
    /* Wait! */
    MPI_Barrier(Universe);
#endif

    /* If we which to precompute the force grid, do it now */
    if (Options.rigidgrid) {
        build3DForceGrid();
    }

    /* Open all the file streams (one for every z point) [ONLY ON ROOT PROCESSOR] */
    if (Me == RootProc) {
        FStreams = (FILE **)malloc((Npoints.z+1)*sizeof(FILE *));
        for (i=0; i<=Npoints.z; ++i) {
            z = Options.zhigh - i*Options.dz;
            if ( Options.gzip == TRUE ) {
                sprintf(outfile,"gzip -6 > scan-%06.3f.dat.gz",z);
                FStreams[i] = popen(outfile,"w");
            }
            else {
                sprintf(outfile,"scan-%06.3f.dat",z);
                FStreams[i] = fopen(outfile,"w");
            }
        }
    }

    /* Note the time */
    gettimeofday(&TimeStart, NULL);

    /* Go home */
    return;
}

/* Now move the tip! */
void moveTip(void) {

    int n, nmax, check;
    int i, ix, iy, iz;
    double x, y, z;
    double angle, minangle, maxforce, dd;
    double e, ediff, eold, fnorm;
    VECTOR f, d;
    VECTOR *ftip;
    int nxy, onproc, bufsize, nsr;
    BUFFER *sendbuf, *recvbuf;
    double checkperc, curperc;

    FILE *tfp;
    char dump[NAME_LENGTH];

    /* Create a force vector (for real time analysis) */
    ftip = (VECTOR *)malloc((Npoints.z+1)*sizeof(VECTOR));

    /* DEBUG DEBUG DEBUG DEBUG DEBUG */
    //Npoints.x = Npoints.y = 10;

    /* Some initialization */
    Ntotal = 0;
    checkperc = curperc = 0.10;
    debugline(RootProc,"Simulation run started now");

    /* Initialize the storage send and receive buffers */
    nxy = (Npoints.x + 1)*(Npoints.y + 1);
    bufsize = Options.bufsize * (Npoints.z+1);
    sendbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    recvbuf = (BUFFER *)malloc(bufsize*sizeof(BUFFER));
    nsr = 0;

    /* Loop x */
    for (ix=0; ix<=Npoints.x; ++ix) {
        x = ix*Options.dx; /* Current x */

        /* Loop y */
        for (iy=0; iy<=Npoints.y; ++iy) {
            y = iy*Options.dy; /* Current y */

            /* Check the progress and report every so often */
            n = ix*(Npoints.y+1) + iy;
            if ( (Me == RootProc) && ((((double)n)/(nxy)) >= curperc) ) {
                debugline(RootProc,"Finished approximately %4.1f %% of the simulation",100*curperc);
                curperc += checkperc;
            }

            /* Compute on which processor this x,y combination should be run */
            onproc = n % NProcessors;
            if (onproc != Me) {
                continue;
            }
            PointsOnProc[Me]++;

            /* Position the tip far above the surface */
            Dummy_pos.x = x;
            Dummy_pos.y = y;
            Dummy_pos.z = Options.zhigh + Options.dz; /* Plus dz to allow for initial subtraction */
            Tip_pos.x = Dummy_pos.x;
            Tip_pos.y = Dummy_pos.y;
            Tip_pos.z = Dummy_pos.z - DummyParams.rmin;

            /* If the molecule is flexible, reset it to its original position, before beginning the approach */
            if (Options.flexible) {
                for (i=0; i<Natoms; ++i) {
                    Surf_pos[i].x = Surf_pos_org[i].x;
                    Surf_pos[i].y = Surf_pos_org[i].y;
                    Surf_pos[i].z = Surf_pos_org[i].z;
                }
            }

            /* Approach and optimize */
            nmax = 0;
            minangle = 9e99;
            maxforce = -9e99;
            for (iz=0; iz<=Npoints.z; ++iz) {
                z = Options.zhigh - iz*Options.dz; /* Current z */

                /* Move tip and dummy atom toward the surface */
                Dummy_pos.z -= Options.dz;
                Tip_pos.z -= Options.dz;

                /* Collect the force */
                ftip[iz] = NULL_vector;

                /* Relax/Minimize the configuration */
                ediff = 5*Options.etol;
                for (n=0; n<Options.maxsteps; ++n) {
                    nmax++;

                    /* Compute all interaction energies */
                    interactTipSurface();
                    interactTipDummy();
                    interactTipHarmonic();

                    /* Total energy and force computed */
                    e = TipSurf_energy + TipDummy_energy + TipHarmonic_energy;
                    f.x = TipSurf_force.x + TipDummy_force.x + TipHarmonic_force.x;
                    f.y = TipSurf_force.y + TipDummy_force.y + TipHarmonic_force.y;
                    f.z = TipSurf_force.z + TipDummy_force.z + TipHarmonic_force.z;
                    fnorm = sqrt(f.x*f.x + f.y*f.y + f.z*f.z);

                    /* Energy difference */
                    if (n>0) {
                        ediff = e - eold;
                    }
                    eold = e;

                    /* Are the forces/energies tolerable */
                    if (Options.minterm == MIN_E) {
                        check = (fabs(ediff)<Options.etol);
                    }
                    else if (Options.minterm == MIN_F) {
                        check = (fabs(fnorm)<Options.ftol);
                    }
                    else if (Options.minterm == MIN_EF) {
                        check = ((fabs(ediff)<Options.etol)&&(fabs(fnorm)<Options.ftol));
                    }
                    if (check) {
                        break;
                    }

                    /* If they are not tolerable, update position of the tip atom based on the force */
                    Tip_pos.x += Options.cfac * f.x;
                    Tip_pos.y += Options.cfac * f.y;
                    Tip_pos.z += Options.cfac * f.z;

                    /* If the molecule is supposed to behave flexible, we need to update its positions too */
                    if (Options.flexible) {
                        updateFlexibleMolecule();
                    }

                } /* End minimization loop */

                /* Compute some other interesting data */
                d.x = Tip_pos.x - x;
                d.y = Tip_pos.y - y;
                d.z = Tip_pos.z - z;
                dd = sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
                angle = atan2(sqrt(d.x*d.x + d.y*d.y),d.z)*(180.0/PI);
                if (angle<minangle) {
                    minangle = angle;
                }
                if (fnorm>maxforce) {
                    maxforce = fnorm;
                }

                /* Collect the final forces on the tip because of the surface */
                ftip[iz].x = TipSurf_force.x;
                ftip[iz].y = TipSurf_force.y;
                ftip[iz].z = TipSurf_force.z;

                /* Store data in send buffers */
                sendbuf[nsr].ix = ix;
                sendbuf[nsr].iy = iy;
                sendbuf[nsr].iz = iz;
                sendbuf[nsr].n  = n;
                sendbuf[nsr].pos.x = x;
                sendbuf[nsr].pos.y = y;
                sendbuf[nsr].pos.z = z;
                sendbuf[nsr].f.x = TipSurf_force.x;
                sendbuf[nsr].f.y = TipSurf_force.y;
                sendbuf[nsr].f.z = TipSurf_force.z;
                sendbuf[nsr].d.x = d.x;
                sendbuf[nsr].d.y = d.y;
                sendbuf[nsr].d.z = d.z;
                sendbuf[nsr].dd = dd;
                sendbuf[nsr].e = TipSurf_energy;
                sendbuf[nsr].angle = angle;
                nsr++;

            } /* End loop in z */

            /* Dump to file (only when the buffer is full, this can and should only happen after a full z approach) */
            if (nsr==bufsize) {
                dumpToFiles(sendbuf,recvbuf,bufsize);
                nsr = 0;
            }

            /* Compute the frequency shift for the given F(z) */
            //computeDeltaF(x,y,ftip);

            if ((ix==(Npoints.x/2)) && (iy==(Npoints.y/2))) {
                sprintf(dump,"test-%d-%d.xyz",ix,iy);
                tfp = fopen(dump,"w");
                fprintf(tfp,"%d\n\n",Natoms);
                for (i=0; i<Natoms; ++i) {
                    fprintf(tfp,"%s %8.4f %8.4f %8.4f\n",Surf_type[i],Surf_pos[i].x,Surf_pos[i].y,Surf_pos[i].z);
                }
                fclose(tfp);
            }

            /* Keep track of counting */
            Ntotal += nmax;

        } /* End loop in y */

    } /* End loop in x */

    /* Dump to file (if it happened that the buffer contains anything at all */
#if !SERIAL
    MPI_Allreduce(&nsr,&n,1,MPI_INT,MPI_SUM,Universe);
#else
    n = nsr;
#endif
    if (n>0) {
        dumpToFiles(sendbuf,recvbuf,nsr);
    }

    /* Say one last thing */
    debugline(RootProc,"Finished %5.1f %% of the simulation",100*curperc);

    /* Go home */
    return;
}

/* Close all the file streams */
void closeUniverse(void) {

    int i;

#if !SERIAL
    /* Wait! */
    MPI_Barrier(Universe);
#endif

    /* Close each separate file stream [ONLY ON ROOT PROCESSOR] */
    if (Me == RootProc) {
        for (i=0; i<=Npoints.z; ++i) {
            if ( Options.gzip == TRUE ) {
                pclose(FStreams[i]);
            }
            else {
                fclose(FStreams[i]);
            }
        }
    }

    /* Go home */
    return;
}

/********************
 ** FINAL THOUGHTS **
 ********************/

void finalize(void) {

    int n, nsum;
    double dtime, timesum;

    /* Note the time */
    gettimeofday(&TimeEnd, NULL);

    /* Time difference */
    dtime = (TimeEnd.tv_sec - TimeStart.tv_sec + (TimeEnd.tv_usec - TimeStart.tv_usec)/1e6);
    timesum = 0.0;
#if !SERIAL
    MPI_Allreduce(&dtime,&timesum,1,MPI_DOUBLE,MPI_SUM,Universe);
#else
    timesum += dtime;
#endif

    /* Collect number of steps from all processors */
    nsum = 0;
#if !SERIAL
    MPI_Allreduce(&Ntotal,&nsum,1,MPI_INT,MPI_SUM,Universe);
#else
    nsum += Ntotal;
#endif

    /* Print some miscelleneous information */
    debugline(RootProc,"Simulation run finished");
    debugline(RootProc,"Statistics:");
    n = (Npoints.x + 1) * (Npoints.y + 1) * (Npoints.z + 1);
    debugline(RootProc,"    Computed %ld tip positions",n);
    debugline(RootProc,"    Needed %ld minimization steps in total",nsum);
    debugline(RootProc,"    Which means approximately %.2f minimization steps per tip position",((double)nsum/n));
    debugline(RootProc,"    The simulation wall time is %.2f seconds",timesum);
    debugline(RootProc,"    The entire simulation took %.2f seconds",dtime);
    debugline(RootProc,"");

    /* Go home */
    return;
}

/***************************
 ** THE PARALLEL UNIVERSE **
 ***************************/

/* Initialize our parallel world */
void openParallelUniverse(int argc, char *argv[]) {

    int i;

#if !SERIAL
    /* Start MPI */
    MPI_Init(&argc,&argv);
#endif

    /* Determine the size of the universe and which processor we are on */
    RootProc = 0;
#if !SERIAL
    Universe = MPI_COMM_WORLD;
    MPI_Comm_rank(Universe,&Me);
    MPI_Comm_size(Universe,&NProcessors);
#else
    Me = 0;
    NProcessors = 1;
#endif

    /* Initialize the checker on how many x,y points for each processor */
    PointsOnProc = (int *)malloc(NProcessors*sizeof(int));
    for (i=0; i<NProcessors; ++i) {
        PointsOnProc[i] = 0;
    }

    /* Go home */
    return;
}

/* Terminate our parallel worlds */
void closeParallelUniverse(void) {

    int i;
    int *pop;

    /* How many x,y points on each processor */
#if !SERIAL
    MPI_Barrier(Universe);
#endif
    pop = (int *)malloc(NProcessors*sizeof(int));
    for (i=0; i<NProcessors; ++i) {
        pop[i] = 0;
    }
#if !SERIAL
    MPI_Allreduce(PointsOnProc,pop,NProcessors,MPI_INT,MPI_SUM,Universe);
#else
    pop[Me] += PointsOnProc[Me];
#endif
    debugline(RootProc,"How many x,y points did each processor handle:");
    for (i=0; i<NProcessors; ++i) {
        debugline(RootProc,"    Processor %2d: %6d x,y points",i,pop[i]);
    }

#if !SERIAL
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

    /* Set up the parallel routines */
    openParallelUniverse(argc,argv);

    /* Initialize the simulation */
    parseCommandLine(argc,argv);    /* Read the command line */
    readInputFile();                            /* Read input file */
    readXYZFile();                              /* Read the XYZ file */
    readParameterFile();                    /* Read the parameter file */

    /* If the molecule is flexible, build the topology */
    if (Options.flexible) {
        buildTopology();
    }

    /* The simulation itself */
    openUniverse();
    moveTip();
    closeUniverse();

    /* Some final thoughts */
    finalize();

    /* And stop the parallel routines properly */
    closeParallelUniverse();

    /* Done */
    return 0;
}
