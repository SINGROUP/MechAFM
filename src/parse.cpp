#include "parse.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.hpp"
#include "messages.hpp"
#include "utility.hpp"
#include "physics.hpp"
#include "simulation.hpp"
#include "vectors.hpp"

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

/* Retrieve the index of the specific atom type */
// int type2num(char *atom) {
    // int i;
    // for (i = 0; i < Natoms; ++i) {
        // if (strcmp(SurfType2Num[i],atom)==0) {
            // break;
        // }
    // }
    // return i;
// }

/* Read stuff from the command line */
void parseCommandLine(int argc, char *argv[], Simulation& simulation) {
    if (simulation.onRootProcessor()) {
        fprintf(stdout,"+ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +\n");
        fprintf(stdout,"|                   Mechanical AFM Model                      |\n");
        fprintf(stdout,"|  Based on: P. Hapala et al, Phys. Rev. B, 90:085421 (2014)  |\n");
        fprintf(stdout,"|            This implemenation by Peter Spijker              |\n");
        fprintf(stdout,"|          2014-2015 (c) Aalto University, Finland            |\n");
        fprintf(stdout,"+ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +\n");
    }
    if ((argc < 2) || (argc > 2)) {
        error(simulation, "Specify an input file to be read!");
    }
    else {
        sprintf(simulation.input_file_name_, "%s", argv[1]);
    }
    return;
}

/* A function to read an input file */
void readInputFile(Simulation& simulation) {

    InputOptions& options = simulation.options_;

    FILE *fp;
    int check;
    char keyword[NAME_LENGTH];
    char value[NAME_LENGTH];
    char line[LINE_LENGTH];
    char tmp_coulomb[NAME_LENGTH], tmp_minterm[NAME_LENGTH];
    char tmp_gzip[NAME_LENGTH], tmp_units[NAME_LENGTH];
    char tmp_flexible[NAME_LENGTH], tmp_rigidgrid[NAME_LENGTH];

    /* Initialize the mandatory options */
    sprintf(options.xyzfile, "");
    sprintf(options.paramfile, "");
    sprintf(options.tipatom, "");
    sprintf(options.dummyatom, "");
    options.minterm = NOT_SET;

    /* Initialize the other options */
    sprintf(options.planeatom, "");
    options.units = U_KCAL;
    sprintf(tmp_units, "%s" ,"kcal/mol");
    options.coulomb = FALSE;
    options.dx = 0.1;
    options.dy = 0.1;
    options.dz = 0.1;
    options.zlow = 6.0;
    options.zhigh = 10.0;
    options.zplane = NEGVAL;
    options.etol = 0.01;
    options.ftol = 0.01;
    options.cfac = 0.001;
    options.maxsteps = 5000;
    options.bufsize = 1000;
    options.gzip = TRUE;
    options.flexible = FALSE;
    options.rigidgrid = FALSE;

    /* Check if the file exists */
    fp = fopen(simulation.input_file_name_, "r");
    if (fp==NULL) {
        error(simulation, "The file %s does not exist!", simulation.input_file_name_);
    }

    /* Scan the file line by line */
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Get the keyword and convert to uppercase */
        sscanf(line, "%s %s", keyword, value);
        strlow(keyword);
        if (strcmp(keyword, "xyzfile") == 0) {
            sprintf(options.xyzfile, "%s", value);
        }
        else if (strcmp(keyword, "paramfile") == 0) {
            sprintf(options.paramfile, "%s", value);
        }
        else if (strcmp(keyword,"tipatom") == 0) {
            sprintf(options.tipatom, "%s", value);
        }
        else if (strcmp(keyword,"dummyatom") == 0) {
            sprintf(options.dummyatom, "%s", value);
        }
        else if (strcmp(keyword, "planeatom") == 0) {
            sprintf(options.planeatom, "%s", value);
        }
        else if (strcmp(keyword, "dx") == 0) {
            options.dx = atof(value);
        }
        else if (strcmp(keyword, "dy") == 0) {
            options.dy = atof(value);
        }
        else if (strcmp(keyword, "dz") == 0) {
            options.dz = atof(value);
        }
        else if (strcmp(keyword, "zlow") == 0) {
            options.zlow = atof(value);
        }
        else if (strcmp(keyword, "zhigh") == 0) {
            options.zhigh = atof(value);
        }
        else if (strcmp(keyword, "zplane") == 0) {
            options.zplane = atof(value);
        }
        else if (strcmp(keyword, "etol") == 0) {
            options.etol = atof(value);
        }
        else if (strcmp(keyword, "ftol") == 0) {
            options.ftol = atof(value);
        }
        else if (strcmp(keyword, "cfac") == 0) {
            options.cfac = atof(value);
        }
        else if (strcmp(keyword, "maxsteps") == 0) {
            options.maxsteps = atoi(value);
        }
        else if (strcmp(keyword, "bufsize") == 0) {
            options.bufsize = atoi(value);
        }
        else if (strcmp(keyword, "gzip") == 0) {
            if (strcmp(value, "on")==0) {
                options.gzip = TRUE;
            }
            else if (strcmp(value, "off") == 0) {
                options.gzip = FALSE;
            }
            else {
                error(simulation, "Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword, "flexible") == 0) {
            if (strcmp(value, "on") == 0) {
                options.flexible = true;
            }
            else if (strcmp(value, "off") == 0) {
                options.flexible = false;
            }
            else {
                error(simulation, "Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword, "rigidgrid") == 0) {
            if (strcmp(value, "on") == 0) {
                options.rigidgrid = true;
            }
            else if (strcmp(value, "off") == 0) {
                options.rigidgrid = false;
            }
            else {
                error(simulation, "Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword, "coulomb") == 0) {
            if (strcmp(value, "on") == 0) {
                options.coulomb = true;
            }
            else if (strcmp(value, "off") == 0) {
                options.coulomb = false;
            }
            else {
                error(simulation, "Option %s must be either on or off!", keyword);
            }
        }
        else if (strcmp(keyword, "minterm") == 0) {
            if (strcmp(value, "e") == 0) {
                options.minterm = MIN_E;
            }
            else if (strcmp(value, "f") == 0) {
                options.minterm = MIN_F;
            }
            else if (strcmp(value, "ef") == 0) {
                options.minterm = MIN_EF;
            }
            else {
                error(simulation, "Option %s must be either e, f or ef!", keyword);
            }
            sprintf(tmp_minterm, "%s", value);
        }
        else if (strcmp(keyword, "units") == 0) {
            if (strcmp(value, "kcal/mol") == 0) {
                options.units = U_KCAL;
            }
            else if (strcmp(value, "kcal") == 0) {
                options.units = U_KCAL;
            }
            else if (strcmp(value, "kJ/mol") == 0) {
                options.units = U_KJ;
            }
            else if (strcmp(value, "kJ") == 0) {
                options.units = U_KJ;
            }
            else if (strcmp(value, "eV") == 0) {
                options.units = U_EV;
            }
            else {
                error(simulation, "Option %s must be either kcal/mol (default), kJ/mol or eV!", keyword);
            }
            sprintf(tmp_units, "%s", value);
        }
        else {
            error(simulation, "Unknown option %s!", keyword);
        }
    }

    /* Check if all necessary options are initialized */
    if (strcmp(options.xyzfile, "") == 0) {
        error(simulation, "Specify at least an xyzfile!");
    }
    if (strcmp(options.paramfile, "") == 0) {
        error(simulation, "Specify at least a parameter file!");
    }
    if (strcmp(options.tipatom, "") == 0) {
        error(simulation, "Specify at least a tip atom!");
    }
    if (strcmp(options.dummyatom, "") == 0) {
        error(simulation, "Specify at least a dummy atom!");
    }
    if (options.minterm == NOT_SET) {
        error(simulation, "Specify at least a minimization termination criterion (e, f, or ef)!");
    }
    if ((strcmp(options.planeatom, "") == 0) && (options.zplane <= NEGVAL)) {
        error(simulation, "Specify at least a plane or a plane atom!");
    }
    else if ((strcmp(options.planeatom, "") != 0) && (options.zplane > NEGVAL)) {
        error(simulation, "Specify only a plane or a plane atom!");
    }

    /* Close file */
    fclose(fp);

    /* Set some useful thingies */
    if (options.coulomb) {
        sprintf(tmp_coulomb, "%s", "on");
    }
    else {
        sprintf(tmp_coulomb, "%s", "off");
    }
    if (options.gzip) {
        sprintf(tmp_gzip, "%s", "on");
    }
    else {
        sprintf(tmp_gzip, "%s", "off");
    }
    if (options.flexible) {
        sprintf(tmp_flexible, "%s", "on");
    }
    else {
        sprintf(tmp_flexible, "%s", "off");
    }
    if (options.rigidgrid) {
        sprintf(tmp_rigidgrid, "%s", "on");
    }
    else {
        sprintf(tmp_rigidgrid, "%s", "off");
    }

    /* Do some sanity checking */
    if ((options.rigidgrid) && (options.flexible)) {
        error(simulation, "Cannot use a flexible molecule with a static force grid!");
    }

    // [> Set function pointers <]
    // if (options.rigidgrid) {
        // interactTipSurface = interactTipSurfaceFromGrid;
    // }
    // else {
        // interactTipSurface = interactTipSurfaceDirectly;
    // }

    /* Talk to me */
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "Input settings for %s:", simulation.input_file_name_);
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "xyzfile:           %-s", options.xyzfile);
    debugline(simulation, simulation.root_processor_, "paramfile:         %-s", options.paramfile);
    debugline(simulation, simulation.root_processor_, "tipatom:           %-s", options.tipatom);
    debugline(simulation, simulation.root_processor_, "dummyatom:         %-s", options.dummyatom);
    if (strcmp(options.planeatom, "")!=0) {
        debugline(simulation, simulation.root_processor_, "planeatom:         %-s", options.planeatom);
    }
    else if (options.zplane > NEGVAL) {
        debugline(simulation, simulation.root_processor_, "zplane:            %-8.4f", options.zplane);
    }
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "units:             %-s", tmp_units);
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "minterm:           %-s", tmp_minterm);
    debugline(simulation, simulation.root_processor_, "etol:              %-8.4f", options.etol);
    debugline(simulation, simulation.root_processor_, "ftol:              %-8.4f", options.ftol);
    debugline(simulation, simulation.root_processor_, "cfac:              %-8.4f", options.cfac);
    debugline(simulation, simulation.root_processor_, "maxsteps:          %-8d", options.maxsteps);
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "zhigh:             %-8.4f", options.zhigh);
    debugline(simulation, simulation.root_processor_, "zlow:              %-8.4f", options.zlow);
    debugline(simulation, simulation.root_processor_, "dx:                %-8.4f", options.dx);
    debugline(simulation, simulation.root_processor_, "dy:                %-8.4f", options.dy);
    debugline(simulation, simulation.root_processor_, "dz:                %-8.4f", options.dz);
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "coulomb:           %-s", tmp_coulomb);
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "flexible:          %-s", tmp_flexible);
    debugline(simulation, simulation.root_processor_, "rigidgrid:         %-s", tmp_rigidgrid);
    debugline(simulation, simulation.root_processor_, "");
    debugline(simulation, simulation.root_processor_, "bufsize:           %-8d", options.bufsize);
    debugline(simulation, simulation.root_processor_, "gzip:              %-s", tmp_gzip);
    debugline(simulation, simulation.root_processor_, "");
    return;
}

/* Read the XYZ file */
void readXYZFile(Simulation& simulation) {

    InputOptions& options = simulation.options_;
    System& system = simulation.system;

    FILE *fp;
    int firstline, realxyz;
    char line[LINE_LENGTH], value[NAME_LENGTH];
    double fdump;

    /* Read the file once, to determine the number of atoms */
    fp = fopen(options.xyzfile, "r");
    if (fp == NULL) {
        error(simulation, "No such file: %s!", options.xyzfile);
    }
    system.n_atoms_ = 0;
    firstline = TRUE;
    realxyz = FALSE;
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* If the first line contains an integer, it is most likely a proper XYZ file */
        if (firstline) {
            sscanf(line, "%s", value);
            if (isint(value)) {
                system.n_atoms_ = atoi(value);
                realxyz = TRUE;
                break;
            }
            firstline = FALSE;
        }
        /* Otherwise continue reading the lines and count how many atoms we have */
        if (!firstline) {
            /* Count useful lines */
            system.n_atoms_++;
        }
    }
    rewind(fp);

    //if (Me == RootProc) fprintf(stdout,"\n*** natoms = %d\n\n",Natoms);

    /* Initialize the global surface atoms vector and lists */
    system.positions_.reserve(system.n_atoms_);
    system.charges_.reserve(system.n_atoms_);
    system.masses_.reserve(system.n_atoms_);
    system.types_.reserve(system.n_atoms_);
    system.fixed_.reserve(system.n_atoms_);
    system.velocities_.reserve(system.n_atoms_);
    system.forces_.reserve(system.n_atoms_);
    for (int i = 0; i < system.n_atoms_; ++i) {
        system.masses_.push_back(0);
    }

    /* Read each line and store the data */
    int n = 0, ncols = 0;
    char dump[LINE_LENGTH], *pch;
    firstline = TRUE;
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        /* If it is a real XYZ file, skip the first two lines */
        if ((realxyz) && (n < 2)) {
            n++;
            continue;
        }
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Based on the first line with actual atom information, determine how many columns there are */
        if (firstline) {
            strcpy(dump, line);
            pch = strtok(dump, " \t\n\r\f");
            while (pch != NULL) {
                pch = strtok(NULL, " \t\n\r\f");
                ncols++;
            }
            firstline = FALSE;
        }
        char type[ATOM_LENGTH];
        double x, y, z, q = 0;
        int fixed = 0;
        if (ncols == 4) {
            sscanf(line, "%s %lf %lf %lf", type, &x, &y, &z);
        } else if (ncols == 5) {
            sscanf(line, "%s %lf %lf %lf %lf", type, &x, &y, &z, &q);
        } else if (ncols == 6) {
            sscanf(line, "%s %lf %lf %lf %lf %d", type, &x, &y, &z, &q, &fixed);
        } else {
            continue;
        }
        system.types_.push_back(std::string(type));
        system.positions_.push_back(Vec3d(x, y, z));
        system.charges_.push_back(q);
        system.fixed_.push_back(fixed);
    }

    for (int i = 0; i < system.n_atoms_; ++i) {
        if (system.fixed_[i] > 0) {
            system.n_fixed_++;
        }
    }

    //if ( Me == RootProc ) {
    //  fprintf(stdout,"\n[%d] *** ncols = %d\n\n",Me,ncols);
    //  for (i=0; i<Natoms; ++i) {
    //  fprintf(stdout,"[%d] >>> %s %lf %lf %lf %lf\n", Me, Surf_type[i], Surf_pos[i].x, Surf_pos[i].y, Surf_pos[i].z, Surf_q[i]);
    //  }
    //}

    setSystemZ(simulation);

    return;
}

void setSystemZ(Simulation& simulation) {
    /* Put the plane atoms at 0 (in z) or as specified separately*/
    int nplaneatoms;
    double avgz;
    InputOptions& options = simulation.options_;
    System& system = simulation.system;

    if (strcmp(options.planeatom, "") != 0) {
        avgz = 0.0;
        nplaneatoms = 0;
        for (int i = 0; i < system.n_atoms_; ++i) {
            if (strcmp(system.types_[i].c_str(), options.planeatom)==0) {
                nplaneatoms++;
                avgz += system.positions_[i].z;
            }
        }
        avgz /= nplaneatoms;
        for (int i = 0; i < system.n_atoms_; ++i) {
            system.positions_[i].z -= avgz;
        }
    }
    else if (options.zplane > NEGVAL) {
        avgz = 0.0;
        for (int i = 0; i < system.n_atoms_; ++i) {
            avgz += system.positions_[i].z;
        }
        avgz /= system.n_atoms_;
        for (int i = 0; i < system.n_atoms_; ++i) {
            system.positions_[i].z -= avgz;
            system.positions_[i].z += options.zplane;
        }
    }
}

void centerSystem(Simulation& simulation){
    InputOptions& options = simulation.options_;
    System& system = simulation.system;

    if (strcmp(options.planeatom, "") != 0) {
        double avgx = 0.0;
        double avgy = 0.0;
        int nplaneatoms = 0;
        for (int i = 0; i < system.n_atoms_; ++i) {
            if (strcmp(system.types_[i].c_str(), options.planeatom) == 0) {
                nplaneatoms++;
                avgx += system.positions_[i].x;
                avgy += system.positions_[i].y;
            }
        }
        avgx /= nplaneatoms;
        avgy /= nplaneatoms;
        double dx = (simulation.box_.x / 2) - avgx;
        double dy = (simulation.box_.y / 2) - avgy;
        for (int i = 0; i < system.n_atoms_; ++i) {
            system.positions_[i].x += dx;
            system.positions_[i].y += dy;
        }
    }
}

/* Read the parameter file */
void readParameterFile(Simulation& simulation) {

    InputOptions& options = simulation.options_;
    System& system = simulation.system;

    FILE *fp;
    char atom[ATOM_LENGTH], keyword[NAME_LENGTH], dump[NAME_LENGTH], line[LINE_LENGTH];
    char atom1[ATOM_LENGTH], atom2[ATOM_LENGTH], style[NAME_LENGTH], *pch;
    double eps, sig, eps_cross, sig_cross, mass, es6, es12;
    double sig6, De, a, re;
    double eps_tip, sig_tip, q_tip, qbase, epermv;
    int check, hcheck, natoms, ncols;
    double chargecheck, qdump;

    /* Initialize the universe */
    simulation.box_ = Vec3d(-1.0);

    /* Open the parameter file */
    fp = fopen(options.paramfile,"r");
    if (fp==NULL) {
        error(simulation, "No parameter file %s found!", options.paramfile);
    }

    /* Scan the parameter file for the universe size and for the tip atom definitions */
    check = FALSE;
    system.n_types_ = 0;
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        /* Skip empty and commented lines */
        if (checkForComments(line)) {
            continue;
        }
        /* Read line to determine keyword */
        sscanf(line, "%s", keyword);
        /* Process the separate keywords */
        if (strcmp(keyword, "box") == 0) {
            Vec3d& box = simulation.box_;
            if ((box.x < 0) && (box.y < 0) && (box.z < 0)) {
                sscanf(line, "%s %lf %lf %lf", dump, &(box.x), &(box.y), &(box.z));
            }
            else {
                error(simulation, "Keyword box cannot be defined more than once in parameter file!");
            }
        }
        if (strcmp(keyword, "atom") == 0) {
            sscanf(line, "%s %s", dump, atom);
            if (strcmp(atom, options.tipatom) == 0) {
                if (check == TRUE) {
                    error(simulation, "Parameters for tip atom can only be specified once!");
                }
                sscanf(line, "%s %s %lf %lf %s %lf", dump, dump, &(eps_tip),
                                                     &(sig_tip), dump,&(q_tip));
                check = TRUE;
            }
            system.n_types_++;
        }
    }
    if (!check) {
        error(simulation, "Parameters for tip atom not defined in parameter file!");
    }
    rewind(fp);

    /* Now we know the size of the universe, put the molecule in the center of it */
    centerSystem(simulation);

    /* Set up the interaction list */
    parseInteractions(simulation, fp);

    fclose(fp);
    return;
}

void parseInteractions(Simulation& simulation, FILE* fp){
    // TipSurfParams = (InteractionList *)malloc(Natoms*sizeof(InteractionList));

    // [> Create a list with all the possible surface atom particle types <]
    // Ntypes -= 2;    [> Subtract the tip and dummy atoms <]
    // if (Ntypes<1) {
        // error(simulation, "Either the tip and/or dummy atom is not specified, or there is no molecule defined. Fix it!");
    // }
    // SurfSurfParams = (InteractionList **)malloc(Ntypes*sizeof(InteractionList *));
    // SurfType2Num = (char **)malloc(Ntypes*sizeof(char *));
    // for (i=0; i<Ntypes; ++i) {
        // SurfSurfParams[i] = (InteractionList *)malloc(Ntypes*sizeof(InteractionList));
        // SurfType2Num[i] = (char *)malloc(ATOM_LENGTH*sizeof(char));
    // }
    // for (j=0, k=0; j<Natoms; ++j) {
        // check = FALSE;
        // for (i=0; i<Ntypes; ++i) {
            // if (strcmp(Surf_type[j],SurfType2Num[i])==0) {
                // check = TRUE;
            // }
        // }
        // if (!check) {
            // strcpy(SurfType2Num[k],Surf_type[j]);
            // k++;
        // }
    // }
    // if (k!=Ntypes) {
        // error(simulation, "Lost an atom type somewhere along the way");
    // }

    // [> The constant part of the Coulomb equation (depends on chosen unit system) <]
    // epermv = 1.0;
    // if (options.units == U_KCAL) {
        // qbase = 332.06371;
    // }
    // else if (options.units == U_KJ) {
        // qbase = 1389.354563;
    // }
    // else if (options.units == U_EV) {
        // qbase = 14.39964901;
    // }
    // qbase /= epermv;

    // [> Quickly check whether charges were read from the XYZ file <]
    // [> NOTE: if only zero charges are specified in the XYZ file, this check sort of fails, <]
    // [>       as the RMSQE will be zero in that case. If you want zero charge, make sure      <]
    // [>       the charges in the parameter file are also set to zero! <]
    // chargecheck = 0.0;
    // for (i=0; i<Natoms; ++i) {
        // chargecheck += (Surf_q[i]*Surf_q[i]);
    // }
    // if (fabs(chargecheck)<TOLERANCE) {
        // warning("The RMSQ-error for the charges read from the XYZ-file is zero. Charges will be read from the parameter file. If you want zero charge, set the charge in the parameter file to zero.");
    // }

    // [> Read the parameter file again, but this time, create the interaction list
         // for the surface atoms, for the dummy atom, and also for the harmonic spring */
    // check = hcheck = FALSE;
    // natoms = 0;
    // while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        // [> Skip empty and commented lines <]
        // if (checkForComments(line)) {
            // continue;
        // }
        // [> Read line to determine keyword <]
        // sscanf(line,"%s",keyword);
        // if (strcmp(keyword,"atom")==0) {
            // sscanf(line,"%s %s %lf %lf %lf %lf",dump,atom,&(eps),&(sig),&(mass),&(qdump));
            // [> Loop all atoms in the surface and check if they match this parameter set <]
            // for (i=0; i<Natoms; ++i) {
                // if (strcmp(Surf_type[i],atom)==0) {
                    // natoms++;
                    // eps_cross = mixeps(eps,eps_tip);
                    // sig_cross = mixsig(sig,sig_tip);                            [> To power 1 <]
                    // sig_cross = (sig_cross*sig_cross*sig_cross);    [> To power 3 <]
                    // sig_cross *= sig_cross;                                             [> To power 6 <]
                    // TipSurfParams[i].es12 = 4 * eps_cross * sig_cross * sig_cross;
                    // TipSurfParams[i].es6    = 4 * eps_cross * sig_cross;
                    // if (fabs(chargecheck)<TOLERANCE) {
                        // Surf_q[i] = qdump;
                    // }
                    // TipSurfParams[i].qq     = qbase * q_tip * Surf_q[i];
                    // TipSurfParams[i].morse = FALSE;
                    // Surf_mass[i] = mass;
                    // [> While we read these parameters and assign the tip-molecule interactions,
                        // we can also build the molecule-molecule interactions (at least the diagonal) */
                    // k = type2num(atom);
                    // SurfSurfParams[k][k].eps = eps;
                    // SurfSurfParams[k][k].sig = sig;
                // }
            // }
            // [> We found a dummy atom in the parameter list <]
            // if (strcmp(options.dummyatom,atom)==0) {
                // if (check == TRUE) {
                    // error(simulation, "Parameters for dummy atom can only be specified once!");
                // }
                // check = TRUE;
                // eps_cross = mixeps(eps,eps_tip);
                // sig_cross = mixsig(sig,sig_tip);                            [> To power 1 <]
                // sig_cross = (sig_cross*sig_cross*sig_cross);    [> To power 3 <]
                // sig_cross *= sig_cross;                                             [> To power 6 <]
                // DummyParams.es12 = 4 * eps_cross * sig_cross * sig_cross;
                // DummyParams.es6  = 4 * eps_cross * sig_cross;
                // DummyParams.qq   = 0.0; [> Ignore Coulomb interaction between tip and dummy <]
                // DummyParams.rmin = mixsig(sig,sig_tip) * SIXTHRT2; [> Needed for tip positioning <]
                // DummyParams.morse = FALSE;
            // }
        // }
        // if (strcmp(keyword,"harm")==0) {
            // if (hcheck == TRUE) {
                // error(simulation, "Parameters for harmonic spring can only be specified once!");
            // }
            // sscanf(line,"%s %s %lf %lf",dump,atom,&(Harmonic.k),&(Harmonic.r0));
            // Harmonic.morse = FALSE;
            // if (strcmp(atom,options.tipatom)!=0) {
                // error(simulation, "Harmonic spring should be defined on tip atom!");
            // }
            // hcheck = TRUE;
        // }
    // }
    // if (natoms != Natoms) {
        // error(simulation, "Not all atoms have been assigned parameters! (%d/%d)",natoms,Natoms);
    // }
    // if (check == FALSE) {
        // error(simulation, "Parameters for dummy atom not defined in parameter file!");
    // }
    // if (hcheck == FALSE) {
        // error(simulation, "No harmonic spring parameters found in parameter file!");
    // }

    // [> Build the entire interaction matrix for surface-surface interactions <]
    // for (i=0; i<Ntypes; ++i) {
        // for (j=0; j<Ntypes; ++j) {
            // eps_cross = mixeps(SurfSurfParams[i][i].eps,SurfSurfParams[j][j].eps);
            // sig_cross = mixsig(SurfSurfParams[i][i].sig,SurfSurfParams[j][j].sig);
            // // sig_cross = mixsig(sig,sig_tip);                            [> To power 1 <]
            // sig_cross = (sig_cross*sig_cross*sig_cross);    [> To power 3 <]
            // sig_cross *= sig_cross;                                             [> To power 6 <]
            // SurfSurfParams[i][j].es12 = 4 * eps_cross * sig_cross * sig_cross;
            // SurfSurfParams[i][j].es6    = 4 * eps_cross * sig_cross;
            // SurfSurfParams[j][i].es12 = SurfSurfParams[i][j].es12;
            // SurfSurfParams[j][i].es6    = SurfSurfParams[i][j].es6;
            // SurfSurfParams[i][j].morse = FALSE;
            // SurfSurfParams[j][i].morse = SurfSurfParams[i][j].morse;
        // }
    // }

    // [> Finally check for pair style overwrites <]
    // rewind(fp);
    // while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        // [> Skip empty and commented lines <]
        // if (checkForComments(line)) {
            // continue;
        // }
        // [> Read line to determine keyword <]
        // sscanf(line,"%s",keyword);
        // if (strcmp(keyword,"pair_ovwrt")==0) {
            // sscanf(line,"%s %s %s %s",dump,atom1,atom2,style);
            // [> Check number of columns <]
            // ncols = 0;
            // strcpy(dump,line);
            // pch = strtok(dump," \t\n\r\f");
            // while (pch!=NULL) {
                // pch = strtok(NULL," \t\n\r\f");
                // ncols++;
            // }
            // [> Lennard-Jones potential <]
            // if (strcmp(style,"lj")==0) {
                // [> Check if number of columns is correct for LJ <]
                // if (ncols>(4+2)) {
                    // error(simulation, "Only two parameters (eps,sig) allowed for LJ");
                // }
                // [> Read overwrite parameters <]
                // sscanf(line,"%s %s %s %s %lf %lf",dump,dump,dump,dump,&(eps),&(sig));
                // sig6 = sig*sig*sig;
                // sig6 *= sig6;
                // es12 = 4 * eps * sig6 * sig6;
                // es6  = 4 * eps * sig6;
            // }
            // [> Morse potential <]
            // else if (strcmp(style,"morse")==0) {
                // [> Check if number of columns is correct for LJ <]
                // if (ncols>(4+3)) {
                    // error(simulation, "Only three parameters (De,a,re) allowed for Morse");
                // }
                // [> Read overwrite parameters <]
                // sscanf(line,"%s %s %s %s %lf %lf %lf",dump,dump,dump,dump,&(De),&(a),&(re));
            // }
            // [> Try and catch <]
            // else {
                // error(simulation, "Unknown pair style overwrite '%s'",style);
            // }
            // [> Store the changes <]
            // if ((strcmp(atom1,options.tipatom)==0) || (strcmp(atom2,options.tipatom)==0)) {
                // if (strcmp(atom1,options.tipatom)==0) {
                    // strcpy(atom,atom2);
                // }
                // if (strcmp(atom2,options.tipatom)==0) {
                    // strcpy(atom,atom1);
                // }
                // [> Lennard-Jones potential <]
                // if (strcmp(style,"lj")==0) {
                    // for (i=0; i<Natoms; ++i) {
                        // if (strcmp(Surf_type[i],atom)==0) {
                            // TipSurfParams[i].es12 = es12;
                            // TipSurfParams[i].es6    = es6;
                            // TipSurfParams[i].morse = FALSE;
                        // }
                    // }
                // }
                // [> Morse potential <]
                // if (strcmp(style,"morse")==0) {
                    // for (i=0; i<Natoms; ++i) {
                        // if (strcmp(Surf_type[i],atom)==0) {
                            // TipSurfParams[i].De = De;
                            // TipSurfParams[i].a  = a;
                            // TipSurfParams[i].re = re;
                            // TipSurfParams[i].morse = TRUE;
                        // }
                    // }
                // }
            // }
            // else {
                // i = type2num(atom1);
                // j = type2num(atom2);
                // if (strcmp(style,"lj")==0) {
                    // SurfSurfParams[i][j].es12 = es12;
                    // SurfSurfParams[i][j].es6 = es6;
                    // SurfSurfParams[j][i].es12 = SurfSurfParams[i][j].es12;
                    // SurfSurfParams[j][i].es6 = SurfSurfParams[i][j].es6;
                    // SurfSurfParams[i][j].morse = FALSE;
                    // SurfSurfParams[j][i].morse = SurfSurfParams[i][j].morse;
                // }
                // if (strcmp(style,"morse")==0) {
                    // SurfSurfParams[i][j].De = De;
                    // SurfSurfParams[i][j].a  = a;
                    // SurfSurfParams[i][j].re = re;
                    // SurfSurfParams[j][i].De = SurfSurfParams[i][j].De;
                    // SurfSurfParams[j][i].a  = SurfSurfParams[i][j].a;
                    // SurfSurfParams[j][i].re = SurfSurfParams[i][j].re;
                    // SurfSurfParams[i][j].morse = TRUE;
                    // SurfSurfParams[j][i].morse = SurfSurfParams[i][j].morse;
                // }
            // }
        // }
    // }

    /* Debug printing */
    //for (i=0; i<Ntypes; ++i) {
    //  for (j=0; j<Ntypes; ++j) {
    //      fprintf(stdout,"%d %d - %8.4f %8.4f - %8.4f %8.4f %8.4f\n",i,j,SurfSurfParams[i][j].es12,SurfSurfParams[i][j].es6,SurfSurfParams[i][j].De,SurfSurfParams[i][j].a,SurfSurfParams[i][j].re);
    //  }
    //}
}
