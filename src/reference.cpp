            // [> Loop all atoms in the surface and check if they match this parameter set <]
            for (i=0; i<Natoms; ++i) {
                if (strcmp(Surf_type[i],atom)==0) {
                    natoms++;
                    eps_cross = mixeps(eps,eps_tip);
                    sig_cross = mixsig(sig,sig_tip);                            // [> To power 1 <]
                    sig_cross = (sig_cross*sig_cross*sig_cross);    // [> To power 3 <]
                    sig_cross *= sig_cross;                                             // [> To power 6 <]
                    TipSurfParams[i].es12 = 4 * eps_cross * sig_cross * sig_cross;
                    TipSurfParams[i].es6    = 4 * eps_cross * sig_cross;
                    if (fabs(chargecheck)<TOLERANCE) {
                        Surf_q[i] = qdump;
                    }
                    TipSurfParams[i].qq     = qbase * q_tip * Surf_q[i];
                    TipSurfParams[i].morse = FALSE;
                    Surf_mass[i] = mass;
                    // While we read these parameters and assign the tip-molecule interactions,
                    // we can also build the molecule-molecule interactions (at least the diagonal)
                    k = type2num(atom);
                    SurfSurfParams[k][k].eps = eps;
                    SurfSurfParams[k][k].sig = sig;
                }
            }
            if (strcmp(options.dummyatom,atom)==0) {
                if (check == TRUE) {
                    error(simulation, "Parameters for dummy atom can only be specified once!");
                }
                check = TRUE;
                eps_cross = mixeps(eps,eps_tip);
                sig_cross = mixsig(sig,sig_tip);                            // [> To power 1 <]
                sig_cross = (sig_cross*sig_cross*sig_cross);    // [> To power 3 <]
                sig_cross *= sig_cross;                                             // [> To power 6 <]
                DummyParams.es12 = 4 * eps_cross * sig_cross * sig_cross;
                DummyParams.es6  = 4 * eps_cross * sig_cross;
                DummyParams.qq   = 0.0; // [> Ignore Coulomb interaction between tip and dummy <]
                DummyParams.rmin = mixsig(sig,sig_tip) * SIXTHRT2; // [> Needed for tip positioning <]
                DummyParams.morse = FALSE;
            }
    }


    // Read the parameter file again, but this time, create the interaction list
    // for the surface atoms, for the dummy atom, and also for the harmonic spring
    check = hcheck = FALSE;
    natoms = 0;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        // [> Skip empty and commented lines <]
        if (checkForComments(line)) {
            continue;
        }
        // [> Read line to determine keyword <]
        sscanf(line,"%s",keyword);
        if (strcmp(keyword,"atom")==0) {
        }
    if (!check) {
        error(simulation, "Parameters for tip atom not defined in parameter file!");
    }
    if (natoms != Natoms) {
        error(simulation, "Not all atoms have been assigned parameters! (%d/%d)",natoms,Natoms);
    }
    if (check == FALSE) {
        error(simulation, "Parameters for dummy atom not defined in parameter file!");
    }

    // [> Build the entire interaction matrix for surface-surface interactions <]
    for (i=0; i<Ntypes; ++i) {
        for (j=0; j<Ntypes; ++j) {
            eps_cross = mixeps(SurfSurfParams[i][i].eps,SurfSurfParams[j][j].eps);
            sig_cross = mixsig(SurfSurfParams[i][i].sig,SurfSurfParams[j][j].sig);
            // sig_cross = mixsig(sig,sig_tip);                            // [> To power 1 <]
            sig_cross = (sig_cross*sig_cross*sig_cross);    // [> To power 3 <]
            sig_cross *= sig_cross;                                             // [> To power 6 <]
            SurfSurfParams[i][j].es12 = 4 * eps_cross * sig_cross * sig_cross;
            SurfSurfParams[i][j].es6    = 4 * eps_cross * sig_cross;
            SurfSurfParams[j][i].es12 = SurfSurfParams[i][j].es12;
            SurfSurfParams[j][i].es6    = SurfSurfParams[i][j].es6;
            SurfSurfParams[i][j].morse = FALSE;
            SurfSurfParams[j][i].morse = SurfSurfParams[i][j].morse;
        }
    }
                // [> Read overwrite parameters <]
                sscanf(line, "%s %s %s %s %lf %lf", dump, dump, dump, dump, &(eps),&(sig));
                sig6 = sig*sig*sig;
                sig6 *= sig6;
                es12 = 4 * eps * sig6 * sig6;
                es6  = 4 * eps * sig6;
            // [> Store the changes <]
            if ((strcmp(atom1,options.tipatom)==0) || (strcmp(atom2,options.tipatom)==0)) {
                if (strcmp(atom1,options.tipatom)==0) {
                    strcpy(atom,atom2);
                }
                if (strcmp(atom2,options.tipatom)==0) {
                    strcpy(atom,atom1);
                }
                // [> Lennard-Jones potential <]
                if (strcmp(style,"lj")==0) {
                    for (i=0; i<Natoms; ++i) {
                        if (strcmp(Surf_type[i],atom)==0) {
                            TipSurfParams[i].es12 = es12;
                            TipSurfParams[i].es6    = es6;
                            TipSurfParams[i].morse = FALSE;
                        }
                    }
                }
                // [> Morse potential <]
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

