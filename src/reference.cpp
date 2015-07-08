void dumpToFiles(BUFFER *sendbuf, BUFFER *recvbuf, int bufsize) {

    int i, f, nsr, *curbufsize, *lcbs;
#if MPI_BUILD
    MPI_Status mpistatus;
#endif

    [> Build an array of the current buffer size for broadcast <]
    curbufsize = (int *)malloc(NProcessors*sizeof(int));
    lcbs = (int *)malloc(NProcessors*sizeof(int));
    for (i=0; i<NProcessors; ++i) {
        lcbs[i] = 0;
    }
    lcbs[Me] = bufsize;
#if MPI_BUILD
    MPI_Allreduce(lcbs,curbufsize,NProcessors,MPI_INT,MPI_SUM,Universe);
#else
    curbufsize[Me] = lcbs[Me];
#endif

#if MPI_BUILD
    [> Send the data to the root processor <]
    if (Me != RootProc) {
        MPI_Send(sendbuf,curbufsize[Me]*sizeof(BUFFER),MPI_CHAR,RootProc,0,Universe);
    }
    [> Receive the date from the daughter processors and write to file <]
    else {
#endif
        [> Loop the processors <]
        for (i=0; i<NProcessors; ++i) {
            [> On the main processor we have to copy the data only, no send and receive <]
            if (i==0) {
                recvbuf = sendbuf;
            }
#if MPI_BUILD
            [> For all other processors we need to receive the data <]
            else {
                MPI_Recv(recvbuf,curbufsize[i]*sizeof(BUFFER),MPI_CHAR,i,0,Universe,&mpistatus);
            }
#endif
            [> Write data to file (only the root processor can do this) <]
            [> PLEASE NOTE: DATA IS SENT IN STRIPED FORM, THEY ARE NOT ORDERED! <]
            for (nsr=0; nsr<curbufsize[i]; ++nsr) {
                f = recvbuf[nsr].iz;
                [> The file buffer can be a gzip pipe or an ASCII file stream <]
                fprintf(FStreams[f],"%d %d %d ",recvbuf[nsr].iz,recvbuf[nsr].ix,recvbuf[nsr].iy);
                fprintf(FStreams[f],"%6.3f %6.3f %6.3f ",recvbuf[nsr].pos.x,recvbuf[nsr].pos.y,recvbuf[nsr].pos.z);
                fprintf(FStreams[f],"%8.4f %8.4f %8.4f ",recvbuf[nsr].f.x,recvbuf[nsr].f.y,recvbuf[nsr].f.z);
                fprintf(FStreams[f],"%6.3f %6.3f %6.3f ",recvbuf[nsr].d.x,recvbuf[nsr].d.y,recvbuf[nsr].d.z);
                fprintf(FStreams[f],"%6.3f %8.4f ",recvbuf[nsr].dd,recvbuf[nsr].angle);
                fprintf(FStreams[f],"%8.4f %d\n",recvbuf[nsr].e,recvbuf[nsr].n);
            }
        }
#if MPI_BUILD
    }
#endif

    [> Get rid of the buffer size broadcast arrays <]
    free(curbufsize);
    free(lcbs);

    [> Go home <]
    return;
}

            // [> Loop all atoms in the surface and check if they match this parameter set <]
            for (i=0; i<Natoms; ++i) {
                if (strcmp(Surf_type[i],atom)==0) {
                    natoms++;
                    eps_cross = mixeps(eps,eps_tip);
                    sig_cross = mixsig(sig,sig_tip);                            // [> To power 1 <]
                    sig_cross = (sig_cross*sig_cross*sig_cross);    // [> To power 3 <]
                    sig_cross *= sig_cross;                                             // [> To power 6 <]
                    TipSurfParams[i].es12 = 4 * eps_cross * sig_cross * sig_cross;
                    TipSurfParams[i].es6 = 4 * eps_cross * sig_cross;
                    if (fabs(chargecheck)<TOLERANCE) {
                        Surf_q[i] = qdump;
                    }
                    TipSurfParams[i].qq = qbase * q_tip * Surf_q[i];
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

