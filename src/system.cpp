#include "system.hpp"

#include <cmath>
#include <cstdio>

#include "globals.hpp"
#include "messages.hpp"


void System::initialize(int n_atoms) {

    n_atoms_ = n_atoms + 2;  // Add tip and dummy to the atom count
    positions_.assign(n_atoms_, Vec3d());
    velocities_.assign(n_atoms_, Vec3d());
    forces_.assign(n_atoms_, Vec3d());
    energies_.assign(n_atoms_, 0);
    charges_.assign(n_atoms_, 0);
    masses_.assign(n_atoms_, 1);
    fixed_.assign(n_atoms_, 0);
    types_.assign(n_atoms_, "");

    // Dummy is allways fixed and tip is allways free to move
    fixed_[0] = 1;
    fixed_[1] = 0;
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
    vector<Vec3d> forces(n_atoms_, Vec3d(0));
    vector<double> energies(n_atoms_, 0);
    for (const auto& interaction : *interactions_) {
        if (interaction->isTipSurface()) {
            interaction->eval(positions_, forces, energies);
        }
    }
    tip_force = forces[1];
    tip_energy = energies[1];
}


void System::makeXYZFile(string folder) const {
    char file_name[NAME_LENGTH];
    sprintf(file_name, "%sstate_%.1f-%.1f-%.1f.xyz", folder.c_str(), positions_[0].x, positions_[0].y, positions_[0].z);
    FILE* file = fopen(file_name, "w");
    fprintf(file, "%d\n\n", n_atoms_);
    for (int i = 0; i < n_atoms_; ++i) {
        fprintf(file, "%s %8.4f %8.4f %8.4f\n", types_[i].c_str(), positions_[i].x,
                                                positions_[i].y, positions_[i].z);
    }
    fclose(file);
    printf("+- Wrote %s\n", file_name);
}


void System::setMoleculeZ() {
    double min_z = 10e10;
    for (int i = 2; i < n_atoms_; ++i) {
        // Ignore hydrogen when looking for the lowest atom
        if (positions_[i].z < min_z && types_[i] != "H"){
            min_z = positions_[i].z;
        }
    }
    // Put the lowest atom at zero (potential minimum)
    for (int i = 2; i < n_atoms_; ++i) {
        positions_[i].z -= min_z;
    }
    offset_.z = -min_z;
}


void System::centerMolecule(Vec2d pos) {
    double avgx = 0.0;
    double avgy = 0.0;
    for (int i = 2; i < n_atoms_; ++i) {
        avgx += positions_[i].x;
        avgy += positions_[i].y;
    }
    avgx /= n_atoms_ - 2;
    avgy /= n_atoms_ - 2;
    double dx = pos.x - avgx;
    double dy = pos.y - avgy;
    for (int i = 2; i < n_atoms_; ++i) {
        positions_[i].x += dx;
        positions_[i].y += dy;
    }
    offset_.x = dx;
    offset_.y = dy;
}

