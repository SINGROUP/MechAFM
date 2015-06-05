#include "flexible.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.hpp"
#include "messages.hpp"
#include "parse.hpp"

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
