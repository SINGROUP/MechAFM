#include "system.hpp"

#include <cmath>
#include <cstdio>

#include "messages.hpp"


// Initialize the state vectors and make sure dummy and tip are included
void System::initialize(int n_atoms) {

    n_atoms_ = n_atoms + 2;  // Add tip and dummy to the atom count
    positions_.assign(n_atoms_, Vec3d());
    velocities_.assign(n_atoms_, Vec3d());
    forces_.assign(n_atoms_, Vec3d());
    energies_.assign(n_atoms_, 0);
    charges_.assign(n_atoms_, 0);
    masses_.assign(n_atoms_, 0);
    fixed_.assign(n_atoms_, false);
    types_.assign(n_atoms_, "");

    // Dummy is allways fixed and tip is allways free to move
    fixed_[0] = true;
    fixed_[1] = false;
}

OutputData System::getOutput() const {
    OutputData data;
    data.position = positions_[0];
    data.r_vec = positions_[1] - positions_[0];
    data.r = data.r_vec.len();
    data.angle = atan2(data.r_vec.getXY().len(), data.r_vec.z) * (180.0 / PI);
    evalTipSurfaceForces(data.tip_force, data.tip_energy);
    return data;
}

void System::evalTipSurfaceForces(Vec3d& tip_force, double& tip_energy) const {
    vector<Vec3d> forces(n_atoms_);
    vector<double> energies(n_atoms_);
    for (const auto& interaction : *interactions_) {
        if (interaction->isTipSurface()) {
            interaction->eval(positions_, forces, energies);
        }
    }
    tip_force = forces[1];
    tip_energy = energies[1];
}

void System::makeXYZFile() const {
    char file_name[NAME_LENGTH];
    sprintf(file_name, "state_%.1f-%.1f-%.1f.xyz", positions_[0].x, positions_[0].y, positions_[0].z);
    FILE* file = fopen(file_name, "w");
    fprintf(file, "%d\n\n", n_atoms_);
    for (int i = 0; i < n_atoms_; ++i) {
        fprintf(file, "%s %8.4f %8.4f %8.4f\n", types_[i].c_str(), positions_[i].x,
                                                positions_[i].y, positions_[i].z);
    }
    fclose(file);
    printf("+- Wrote %s\n", file_name);
}
