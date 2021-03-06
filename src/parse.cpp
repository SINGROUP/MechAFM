#include "parse.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
    #include <windows.h>
#endif

#include "globals.hpp"
#include "messages.hpp"
#include "simulation.hpp"
#include "utility.hpp"
#include "vectors.hpp"

// Filter out comment and empty lines in a file
bool checkForComments(char* line) {
    bool moveon = false;
    if (line[0] == '#') {
        moveon = true;
    } else if (line[0] == '%') {
        moveon = true;
    } else if (line[0] == '\n') {
        moveon = true;
    }
    return moveon;
}

// Read stuff from the command line
void parseCommandLine(int argc, char* argv[], Simulation& simulation) {
    if (simulation.rootProcess()) {
        fprintf(stdout,"+ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +\n");
        fprintf(stdout,"|                   Mechanical AFM Model                      |\n");
        fprintf(stdout,"|  Based on: P. Hapala et al, Phys. Rev. B, 90:085421 (2014)  |\n");
        fprintf(stdout,"|                  This implementation by                     |\n");
        fprintf(stdout,"|             Peter Spijker and Olli Keisanen                 |\n");
        fprintf(stdout,"|          2014-2015 (c) Aalto University, Finland            |\n");
        fprintf(stdout,"+ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +\n");
    }
    InputOptions& options = simulation.options_;
    if ((argc < 2) ) {
        error("Specify an input file to be read!");
    } else if (argc > 3) {
        error("Too many command line arguments!");
    }
    options.inputfolder = "";
    options.outputfolder = "";
    if (argc >= 2) {
        options.inputfile = argv[1];
#ifdef _WIN32
        size_t path_split = options.inputfile.rfind("\\");
#else
        size_t path_split = options.inputfile.rfind("/");
#endif
        if (path_split != string::npos) {
            options.inputfolder = options.inputfile.substr(0, path_split + 1);
        }
    }
    if (argc == 3) {
        options.outputfolder = argv[2];
#ifdef _WIN32
        if (options.outputfolder[options.outputfolder.size() - 1] != '\\') {
            options.outputfolder += '\\';
        }
#else
        if (options.outputfolder[options.outputfolder.size() - 1] != '/') {
            options.outputfolder += '/';
        }
#endif
        if (simulation.rootProcess()) {
#ifdef _WIN32
            CreateDirectory(options.outputfolder.c_str(), NULL);
#else
            string dir_cmd = "mkdir -p " + options.outputfolder;
            system(dir_cmd.c_str());
#endif
        }
    }
    return;
}

// A function to read an input file
void readInputFile(Simulation& simulation) {
    InputOptions& options = simulation.options_;

    FILE *fp;
    char keyword[NAME_LENGTH];
    char value[NAME_LENGTH];
    char line[LINE_LENGTH];
    char dump[LINE_LENGTH];
    char tmp_coulomb[NAME_LENGTH], tmp_tip_dummy_coulomb[NAME_LENGTH], tmp_minterm[NAME_LENGTH];
    char tmp_gzip[NAME_LENGTH], tmp_statistics[NAME_LENGTH], tmp_units[NAME_LENGTH];
    char tmp_flexible[NAME_LENGTH], tmp_rigidgrid[NAME_LENGTH], tmp_normal[NAME_LENGTH];
    char tmp_use_external_potential[NAME_LENGTH], tmp_vdw_pbc[NAME_LENGTH];

    // Initialize the mandatory options
    options.xyzfile = "";
    options.paramfile = "";
    options.tipatom = "";
    options.dummyatom = "";
    options.minterm = NOT_SET;

    // Initialize the other options
    options.planeatom = "";
    options.units = U_KCAL;
    sprintf(tmp_units, "%s" ,"kcal/mol");
    options.e_potential_file = "";
    options.coulomb = false;
    options.tip_dummy_coulomb = false;
    options.use_external_potential = false;
    options.area = Vec2d(10);
    options.center = Vec2d(-1);
    options.dx = 0.1;
    options.dy = 0.1;
    options.dz = 0.1;
    options.zlow = 6.0;
    options.zhigh = 10.0;
    options.vdw_pbc = false;
    options.cell_a = Vec3d(0);
    options.cell_b = Vec3d(0);
    options.cell_c = Vec3d(0);
    options.normal = NORMAL_Z;
    sprintf(tmp_normal, "%s" ,"z");
    options.etol = 0.01;
    options.ftol = 0.01;
    options.dt = 0.001;
    options.maxsteps = 5000;
    options.bufsize = 1000;
    options.gzip = true;
    options.statistics = false;
    options.flexible = false;
    options.rigidgrid = false;
    options.minimiser_type = FIRE;
    options.integrator_type = MIDPOINT;

    // Check if the file exists
    fp = fopen(options.inputfile.c_str(), "r");
    if (fp == NULL) {
        error("The file %s does not exist!", options.inputfile.c_str());
    }

    // Scan the file line by line
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        // Skip empty and commented lines
        if (checkForComments(line)) {
            continue;
        }
        // Get the keyword and convert to lowercase
        sscanf(line, "%s %s", keyword, value);
        strlow(keyword);
        if (strcmp(keyword, "xyzfile") == 0) {
            options.xyzfile = options.inputfolder + value;
        } else if (strcmp(keyword, "paramfile") == 0) {
            options.paramfile = options.inputfolder + value;
        } else if (strcmp(keyword, "e_potential_file") == 0) {
            options.e_potential_file = options.inputfolder + value;
        } else if (strcmp(keyword, "tipatom") == 0) {
            options.tipatom = value;
        } else if (strcmp(keyword, "dummyatom") == 0) {
            options.dummyatom = value;
        } else if (strcmp(keyword, "area") == 0) {
            sscanf(line, "%s %lf %lf", dump, &(options.area.x), &(options.area.y));
        } else if (strcmp(keyword, "center") == 0) {
            sscanf(line, "%s %lf %lf", dump, &(options.center.x), &(options.center.y));
        } else if (strcmp(keyword, "zhigh") == 0) {
            options.zhigh = atof(value);
        } else if (strcmp(keyword, "zlow") == 0) {
            options.zlow = atof(value);
        } else if (strcmp(keyword, "dx") == 0) {
            options.dx = atof(value);
        } else if (strcmp(keyword, "dy") == 0) {
            options.dy = atof(value);
        } else if (strcmp(keyword, "dz") == 0) {
            options.dz = atof(value);
        } else if (strcmp(keyword, "cell_a") == 0) {
            sscanf(line, "%s %lf %lf %lf", dump, &(options.cell_a.x), &(options.cell_a.y), &(options.cell_a.z));
        } else if (strcmp(keyword, "cell_b") == 0) {
            sscanf(line, "%s %lf %lf %lf", dump, &(options.cell_b.x), &(options.cell_b.y), &(options.cell_b.z));
        } else if (strcmp(keyword, "cell_c") == 0) {
            sscanf(line, "%s %lf %lf %lf", dump, &(options.cell_c.x), &(options.cell_c.y), &(options.cell_c.z));
        } else if (strcmp(keyword, "etol") == 0) {
            options.etol = atof(value);
        } else if (strcmp(keyword, "ftol") == 0) {
            options.ftol = atof(value);
        } else if (strcmp(keyword, "dt") == 0) {
            options.dt = atof(value);
        } else if (strcmp(keyword, "maxsteps") == 0) {
            options.maxsteps = atoi(value);
        } else if (strcmp(keyword, "bufsize") == 0) {
            options.bufsize = atoi(value);
        } else if (strcmp(keyword, "gzip") == 0) {
            if (strcmp(value, "on") == 0) {
                options.gzip = true;
            } else if (strcmp(value, "off") == 0) {
                options.gzip = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if (strcmp(keyword, "statistics") == 0) {
            if (strcmp(value, "on") == 0) {
                options.statistics = true;
            } else if (strcmp(value, "off") == 0) {
                options.statistics = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if (strcmp(keyword, "flexible") == 0) {
            if (strcmp(value, "on") == 0) {
                options.flexible = true;
            } else if (strcmp(value, "off") == 0) {
                options.flexible = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if (strcmp(keyword, "rigidgrid") == 0) {
            if (strcmp(value, "on") == 0) {
                options.rigidgrid = true;
            } else if (strcmp(value, "off") == 0) {
                options.rigidgrid = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if (strcmp(keyword, "coulomb") == 0) {
            if (strcmp(value, "on") == 0) {
                options.coulomb = true;
            } else if (strcmp(value, "off") == 0) {
                options.coulomb = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if (strcmp(keyword, "tip_dummy_coulomb") == 0) {
            if (strcmp(value, "on") == 0) {
                options.tip_dummy_coulomb = true;
            } else if (strcmp(value, "off") == 0) {
                options.tip_dummy_coulomb = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if (strcmp(keyword, "use_external_potential") == 0) {
            if (strcmp(value, "on") == 0) {
                options.use_external_potential = true;
            } else if (strcmp(value, "off") == 0) {
                options.use_external_potential = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if (strcmp(keyword, "vdw_pbc") == 0) {
            if (strcmp(value, "on") == 0) {
                options.vdw_pbc = true;
            } else if (strcmp(value, "off") == 0) {
                options.vdw_pbc = false;
            } else {
                error("Option %s must be either on or off!", keyword);
            }
        } else if(strcmp(keyword, "surface_normal") == 0) {
            strlow(value);
            if (strcmp(value, "z") == 0) {
                options.normal = NORMAL_Z;
            } else if (strcmp(value, "y") == 0) {
                options.normal = NORMAL_Y;
            } else if (strcmp(value, "x") == 0) {
                options.normal = NORMAL_X;
            } else {
                error("Option %s must be either X, Y or Z!", keyword);
            }
            sprintf(tmp_normal, "%s", value);
        } else if (strcmp(keyword, "minterm") == 0) {
            if (strcmp(value, "e") == 0) {
                options.minterm = MIN_E;
            } else if (strcmp(value, "f") == 0) {
                options.minterm = MIN_F;
            } else if (strcmp(value, "ef") == 0) {
                options.minterm = MIN_EF;
            } else {
                error("Option %s must be either e, f or ef!", keyword);
            }
            sprintf(tmp_minterm, "%s", value);
        } else if (strcmp(keyword, "units") == 0) {
            double epermv = 1.0;
            if (strcmp(value, "kcal/mol") == 0) {
                options.units = U_KCAL;
                simulation.interaction_parameters_.qbase = 332.06371 / epermv;
            } else if (strcmp(value, "kcal") == 0) {
                options.units = U_KCAL;
                simulation.interaction_parameters_.qbase = 332.06371 / epermv;
            } else if (strcmp(value, "kJ/mol") == 0) {
                options.units = U_KJ;
                simulation.interaction_parameters_.qbase = 1389.354563 / epermv;
            } else if (strcmp(value, "kJ") == 0) {
                options.units = U_KJ;
                simulation.interaction_parameters_.qbase = 1389.354563 / epermv;
            } else if (strcmp(value, "eV") == 0) {
                options.units = U_EV;
                simulation.interaction_parameters_.qbase = 14.39964901 / epermv;
            } else {
                error("Option %s must be either kcal/mol (default), kJ/mol or eV!", keyword);
            }
            sprintf(tmp_units, "%s", value);
        } else if (strcmp(keyword, "minimiser") == 0) {
            if (strcmp(value, "SD") == 0) {
                options.minimiser_type = STEEPEST_DESCENT;
            } else if (strcmp(value, "FIRE") == 0) {
                options.minimiser_type = FIRE;
            } else {
                error("Unrecognised minimiser type!");
            }
        } else if (strcmp(keyword, "integrator") == 0) {
            if (strcmp(value, "euler") == 0) {
                options.integrator_type = EULER;
            } else if (strcmp(value, "midpoint") == 0) {
                options.integrator_type = MIDPOINT;
            } else if (strcmp(value, "rk4") == 0) {
                options.integrator_type = RK4;
            } else {
                error("Unrecognised integrator type!");
            }
        } else {
            error("Unknown option %s!", keyword);
        }
    }

    // Check if all necessary options are initialized
    if (options.xyzfile == "") {
        error("Specify at least an xyzfile!");
    }
    if (options.paramfile == "") {
        error("Specify at least a parameter file!");
    }
    if (options.tipatom == "") {
        error("Specify at least a tip atom!");
    }
    if (options.dummyatom == "") {
        error("Specify at least a dummy atom!");
    }
    if (options.minterm == NOT_SET) {
        error("Specify at least a minimization termination criterion (e, f, or ef)!");
    }

    fclose(fp);

    // Set some useful thingies
    if (options.center == Vec2d(-1)){
        options.center = options.area / 2;
    }
    if (options.coulomb) {
        sprintf(tmp_coulomb, "%s", "on");
    } else {
        sprintf(tmp_coulomb, "%s", "off");
    }
    if (options.tip_dummy_coulomb) {
        sprintf(tmp_tip_dummy_coulomb, "%s", "on");
    } else {
        sprintf(tmp_tip_dummy_coulomb, "%s", "off");
    }
    if (options.use_external_potential) {
        sprintf(tmp_use_external_potential, "%s", "on");
    } else {
        sprintf(tmp_use_external_potential, "%s", "off");
    }
    if (options.vdw_pbc) {
        sprintf(tmp_vdw_pbc, "%s", "on");
    } else {
        sprintf(tmp_vdw_pbc, "%s", "off");
    }
    if (options.gzip) {
        sprintf(tmp_gzip, "%s", "on");
    } else {
        sprintf(tmp_gzip, "%s", "off");
    }
    if (options.statistics) {
        sprintf(tmp_statistics, "%s", "on");
    } else {
        sprintf(tmp_statistics, "%s", "off");
    }
    if (options.flexible) {
        sprintf(tmp_flexible, "%s", "on");
    } else {
        sprintf(tmp_flexible, "%s", "off");
    }
    if (options.rigidgrid) {
        sprintf(tmp_rigidgrid, "%s", "on");
    } else {
        sprintf(tmp_rigidgrid, "%s", "off");
    }

    // Do some sanity checking
    if ((options.rigidgrid) && (options.flexible)) {
        error("Cannot use a flexible molecule with a static force grid!");
    }
    if (options.coulomb && options.use_external_potential) {
        error("Cannot use Coulomb interaction and external electrostatic potential at the same time!");
    }
    if (options.use_external_potential && options.flexible) {
        error("External potential can be used only for non-flexible systems!");
    }
    if (options.use_external_potential && (options.e_potential_file == "")) {
        error("If you want to use external electrostatic potential, you must specify a file that contains it!");
    }
    if (options.vdw_pbc && options.coulomb) {
        error("Implementation of Coulomb interaction does not support any periodic boundary conditions! Use periodic external electrostatic potential instead.");
    }
    if (options.vdw_pbc && !options.use_external_potential) {
        if (options.cell_a == Vec3d(0) || options.cell_b == Vec3d(0) || options.cell_c == Vec3d(0))
            error("The unit cell vectors must be given if periodic vdW is used.");
    }

    // Talk to me
    pretty_print("");
    pretty_print("Input settings for        %s:", options.inputfile.c_str());
    pretty_print("");
    pretty_print("xyzfile:                  %-s", options.xyzfile.c_str());
    pretty_print("paramfile:                %-s", options.paramfile.c_str());
    pretty_print("tipatom:                  %-s", options.tipatom.c_str());
    pretty_print("dummyatom:                %-s", options.dummyatom.c_str());
    pretty_print("");
    pretty_print("units:                    %-s", tmp_units);
    pretty_print("");
    pretty_print("minterm:                  %-s", tmp_minterm);
    pretty_print("etol:                     %-8.4f", options.etol);
    pretty_print("ftol:                     %-8.4f", options.ftol);
    pretty_print("dt:                       %-8.4f", options.dt);
    pretty_print("maxsteps:                 %-8d", options.maxsteps);
    pretty_print("");
    pretty_print("area:                     %-8.4f %-8.4f", options.area.x, options.area.y);
    pretty_print("center:                   %-8.4f %-8.4f", options.center.x, options.center.y);
    pretty_print("zhigh:                    %-8.4f", options.zhigh);
    pretty_print("zlow:                     %-8.4f", options.zlow);
    pretty_print("dx:                       %-8.4f", options.dx);
    pretty_print("dy:                       %-8.4f", options.dy);
    pretty_print("dz:                       %-8.4f", options.dz);
    pretty_print("");
    pretty_print("surface_normal:           %-s", tmp_normal);
    pretty_print("vdw_pbc:                  %-s", tmp_vdw_pbc);
    if (options.vdw_pbc && !options.use_external_potential) {
        pretty_print("cell_a:                   %-8.4f %-8.4f %-8.4f", options.cell_a.x, options.cell_a.y, options.cell_a.z);
        pretty_print("cell_b:                   %-8.4f %-8.4f %-8.4f", options.cell_b.x, options.cell_b.y, options.cell_b.z);
        pretty_print("cell_c:                   %-8.4f %-8.4f %-8.4f", options.cell_c.x, options.cell_c.y, options.cell_c.z);
    }
    pretty_print("");
    pretty_print("coulomb:                  %-s", tmp_coulomb);
    pretty_print("tip_dummy_coulomb:        %-s", tmp_tip_dummy_coulomb);
    pretty_print("use_external_potential:   %-s", tmp_use_external_potential);
    if (options.use_external_potential)
        pretty_print("e_potential_file:         %-s", options.e_potential_file.c_str());
    pretty_print("");
    pretty_print("flexible:                 %-s", tmp_flexible);
    pretty_print("rigidgrid:                %-s", tmp_rigidgrid);
    pretty_print("");
    switch (options.minimiser_type) {
        case STEEPEST_DESCENT:
            pretty_print("minimiser:         %-s", "SD");
            break;
        case FIRE:
            pretty_print("minimiser:         %-s", "FIRE");
            break;
    }
    switch (options.integrator_type) {
        case EULER:
            pretty_print("integrator:        %-s", "euler");
            break;
        case MIDPOINT:
            pretty_print("integrator:        %-s", "midpoint");
            break;
        case RK4:
            pretty_print("integrator:        %-s", "rk4");
            break;
    }
    pretty_print("");
    pretty_print("bufsize:           %-8d", options.bufsize);
    pretty_print("gzip:              %-s", tmp_gzip);
    pretty_print("statistics:        %-s", tmp_statistics);
    pretty_print("");
    return;
}

// Read the XYZ file
void readXYZFile(Simulation& simulation) {
    InputOptions& options = simulation.options_;
    System& system = simulation.system;

    FILE *fp;
    bool firstline, realxyz;
    char line[LINE_LENGTH], value[NAME_LENGTH];

    // Read the file once, to determine the number of atoms
    fp = fopen(options.xyzfile.c_str(), "r");
    if (fp == NULL) {
        error("No such file: %s!", options.xyzfile.c_str());
    }
    int n_atoms = 0;
    firstline = true;
    realxyz = false;
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        // Skip empty and commented lines
        if (checkForComments(line)) {
            continue;
        }
        // If the first line contains an integer, it is most likely a proper XYZ file
        if (firstline) {
            sscanf(line, "%s", value);
            if (isint(value)) {
                n_atoms = atoi(value);
                realxyz = true;
                break;
            }
            firstline = false;
        }
        // Otherwise continue reading the lines and count how many atoms we have
        if (!firstline) {
            // Count useful lines
            n_atoms++;
        }
    }
    rewind(fp);

    // Initialize the system and set the tip and dummy types which aren't
    // in the parameter file.
    system.initialize(n_atoms);
    system.types_[0] = options.dummyatom;
    system.types_[1] = options.tipatom;

    // Read each line and store the data
    int atom_i = 2;  // Start indexing from 2 so that we don't overwrite the tip or dummy
    int n = 0, ncols = 0;
    char dump[LINE_LENGTH], *pch;
    firstline = true;
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        // If it is a real XYZ file, skip the first two lines
        if ((realxyz) && (n < 2)) {
            n++;
            continue;
        }
        // Skip empty and commented lines
        if (checkForComments(line)) {
            continue;
        }
        // Based on the first line with actual atom information, determine how many columns there are
        if (firstline) {
            strcpy(dump, line);
            pch = strtok(dump, " \t\n\r\f");
            while (pch != NULL) {
                pch = strtok(NULL, " \t\n\r\f");
                ncols++;
            }
            firstline = false;
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
            error("Invalid number of columns in the xyz file.");
        }
        system.types_[atom_i] = std::string(type);
        system.positions_[atom_i] = Vec3d(x, y, z);
        system.charges_[atom_i] = q;
        // Keep atoms fixed if we're not flexible.
        system.fixed_[atom_i] = options.flexible ? fixed : 1;
        atom_i++;
    }

    // Quickly check whether charges were read from the XYZ file
    // NOTE: if only zero charges are specified in the XYZ file, this check sort of fails,
    //       as the RMSQE will be zero in that case. If you want zero charge, make sure
    //       the charges in the parameter file are also set to zero!
    double charge_check = 0;
    for (int i = 0; i < system.n_atoms_; ++i) {
        charge_check += pow(system.charges_[i], 2);
    }
    if (fabs(charge_check) < TOLERANCE) {
        options.xyz_charges = false;
        warning("The RMSQ-error for the charges read from the XYZ-file is zero.\n"
            "+-     Charges will be read from the parameter file.\n"
            "+-     If you want zero charge, set the charge in the parameter file to zero.\n"
            "+-");
    } else {
        options.xyz_charges = true;
    }
    return;
}

// Read the parameter file
void readParameterFile(Simulation& simulation) {
    InputOptions& options = simulation.options_;
    InteractionParameters& parameters = simulation.interaction_parameters_;

    FILE *fp;
    char atom[ATOM_LENGTH], keyword[NAME_LENGTH], dump[NAME_LENGTH], line[LINE_LENGTH];
    char atom1[ATOM_LENGTH], atom2[ATOM_LENGTH], style[NAME_LENGTH], *pch;
    double eps, sig, mass;
    double De, a, re;
    double qdump;

    // Open the parameter file
    fp = fopen(options.paramfile.c_str(), "r");
    if (fp == NULL) {
        error("No parameter file %s found!", options.paramfile.c_str());
    }

    // Scan the parameter file
    bool hcheck = false;
    while (fgets(line, LINE_LENGTH, fp) != NULL) {
        // Skip empty and commented lines
        if (checkForComments(line)) {
            continue;
        }
        // Read line to determine keyword
        sscanf(line, "%s", keyword);
        // Process the separate keywords
        if (strcmp(keyword, "atom") == 0) {
            sscanf(line, "%s %s %lf %lf %lf %lf", dump, atom, &(eps), &(sig), &(mass), &(qdump));
            if (parameters.atom_parameters.count(atom) != 0){
                warning("Parameters for atom type %s defined multiple times.", atom);
            }
            auto& p = parameters.atom_parameters[atom];
            p.eps = eps;
            p.sig = sig;
            p.q = qdump;
            p.mass = mass;
            // Set tip and dummy charges since they won't be read from the xyz file.
            if (options.dummyatom == atom) {
                simulation.system.charges_[0] = qdump;
            }
            if (options.tipatom == atom) {
                simulation.system.charges_[1] = qdump;
            }
        } else if (strcmp(keyword, "harm") == 0) {
            if (hcheck) {
                error("Parameters for harmonic spring can only be specified once!");
            }
            sscanf(line,"%s %lf %lf", dump, &(parameters.tip_dummy_k),
                                            &(parameters.tip_dummy_r0));
            hcheck = true;
        } else if (strcmp(keyword, "pair_ovwrt") == 0) {
            sscanf(line, "%s %s %s %s", dump, atom1, atom2, style);
            // Check number of columns
            int ncols = 0;
            strcpy(dump, line);
            pch = strtok(dump, " \t\n\r\f");
            while (pch != NULL) {
                pch = strtok(NULL, " \t\n\r\f");
                ncols++;
            }
            // Lennard-Jones potential
            if (strcmp(style, "lj") == 0) {
                // Check if number of columns is correct for LJ
                if (ncols != (4 + 2)) {
                    error("Only two parameters (eps, sig) allowed for LJ");
                }
                // Read overwrite parameters
                sscanf(line, "%s %s %s %s %lf %lf", dump, dump, dump, dump, &(eps), &(sig));
                OverwriteParameters p;
                p.atoms = unordered_multiset<string>{atom1, atom2};
                p.eps = eps;
                p.sig = sig;
                p.morse = false;
                parameters.overwrite_parameters.push_back(p);
            }
            // Morse potential
            else if (strcmp(style, "morse") == 0) {
                // Check if number of columns is correct for Morse
                if (ncols != (4 + 3)) {
                    error("Only three parameters (De, a, re) allowed for Morse");
                }
                // Read overwrite parameters
                sscanf(line, "%s %s %s %s %lf %lf %lf", dump, dump, dump, dump, &(De), &(a), &(re));
                OverwriteParameters p;
                p.atoms = unordered_multiset<string>{atom1, atom2};
                p.de = De;
                p.a = a;
                p.re = re;
                p.morse = true;
                parameters.overwrite_parameters.push_back(p);
            }
            // Try and catch
            else {
                error("Unknown pair style overwrite '%s'",style);
            }
        }
    }
    if (hcheck == false) {
        error("No harmonic spring parameters found in parameter file!");
    }

    if (options.flexible) {
        rewind(fp);
        readFlexibleParameters(simulation, fp);
    }

    fclose(fp);
    return;
}

void readFlexibleParameters(Simulation& simulation, FILE* fp) {
    InteractionParameters& parameters = simulation.interaction_parameters_;

    char line[LINE_LENGTH];
    char keyword[NAME_LENGTH], dump[NAME_LENGTH];
    char atom1[NAME_LENGTH], atom2[NAME_LENGTH];

    // Read the parameter file and parse everything we need
    // for modeling a flexible molecule
    double r0;
    bool bcheck, acheck, dcheck, scheck, tcheck;
    bcheck = acheck = dcheck = scheck = false;
    while (fgets(line, LINE_LENGTH, fp)!=NULL) {
        // Skip empty and commented lines
        if (checkForComments(line)) {
            continue;
        }
        // Read line to determine keyword
        sscanf(line,"%s",keyword);
        // The strength of the harmonic bond from the parameter file
        if (strcmp(keyword, "bond") == 0) {
            if (bcheck == true) {
                warning("Parameters for harmonic bond defined multiple times!");
            }
            sscanf(line,"%s %lf", dump, &(parameters.bond_k));
            bcheck = true;
        }
        // The strength of the harmonic angle from the parameter file
        if (strcmp(keyword, "angle") == 0) {
            if (acheck == true) {
                warning("Parameters for harmonic angle defined multiple times!");
            }
            sscanf(line,"%s %lf", dump, &(parameters.angle_k));
            acheck = true;
        }
        // The strength of the harmonic dihedral from the parameter file
        if (strcmp(keyword, "dihedral") == 0) {
            if (dcheck == true) {
                warning("Parameters for harmonic dihedral defined multiple times!");
            }
            sscanf(line,"%s %lf", dump, &(parameters.dihedral_k));
            dcheck = true;
        }
        // The strength of the substrate support from the parameter file
        if (strcmp(keyword, "substrate") == 0) {
            if (scheck == true) {
                warning("Parameters for substrate support defined multiple times!");
            }
            sscanf(line, "%s %lf %lf %lf %lf %lf", dump, &parameters.substrate_eps,
		   &parameters.substrate_sig, &parameters.substrate_lambda, &parameters.substrate_rc, &parameters.substrate_k);
            scheck = true;
        }
        // Collect the bonds, possibly present in the system
        if (strcmp(keyword, "topobond") == 0) {
            tcheck = false;
            sscanf(line, "%s %s %s %lf", dump, atom1, atom2, &(r0));
            unordered_multiset<string> atoms{atom1, atom2};
            for (const auto& bond : parameters.possible_bonds_) {
                if (bond.atoms == atoms){
                    tcheck = true;
                }
            }
            if (tcheck) {
                warning("The topobond for %s and %s is defined at least twice!", atom1, atom2);
            }
            PossibleBond pb;
            pb.atoms = atoms;
            pb.r0 = r0;
            parameters.possible_bonds_.push_back(pb);
        }
    }
    if (bcheck == false) {
        error("No harmonic bond parameters found in parameter file!");
    }
    if (acheck == false) {
        error("No harmonic angle parameters found in parameter file!");
    }
    if (dcheck == false) {
        error("No harmonic dihedral parameters found in parameter file!");
    }
    if (scheck == false) {
        error("No harmonic substrate support parameters found in parameter file!");
    }
}
