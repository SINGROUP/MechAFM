#include "parse.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.hpp"
#include "messages.hpp"
#include "utility.hpp"
#include "physics.hpp"
#include "simulation.hpp"

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
int type2num(char *atom) {
    int i;
    for (i=0; i<Ntypes; ++i) {
        if (strcmp(SurfType2Num[i],atom)==0) {
            break;
        }
    }
    return i;
}

/* Read stuff from the command line */
void parseCommandLine(int argc, char *argv[], Simulation& simulation) {
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
