#include "system.hpp"

#include <cmath>
#include <cstdio>

#include "globals.hpp"
#include "messages.hpp"


void System::initialize(int n_atoms) {
    tip_pbc_ = false;
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
    data.position.x = real_tip_xy_.x;
    data.position.y = real_tip_xy_.y;
    data.position.z = positions_[0].z;
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
    sprintf(file_name, "%sstate_%.1f-%.1f-%.1f.xyz", folder.c_str(), real_tip_xy_.x, real_tip_xy_.y, positions_[0].z);
    FILE* file = fopen(file_name, "w");
    fprintf(file, "%d\n\n", n_atoms_);
    for (int i = 0; i < n_atoms_; ++i) {
        fprintf(file, "%s %8.4f %8.4f %8.4f\n", types_[i].c_str(), positions_[i].x,
                                                positions_[i].y, positions_[i].z);
    }
    fclose(file);
    printf("+- Wrote %s\n", file_name);
}


void System::rotateCoordAxes(const string& new_coord_sequence) {
    vector<Vec3d> new_positions;
    new_positions.resize(positions_.size());
    
    if (new_coord_sequence == "ZXY") {
        for (unsigned int ia = 0; ia < positions_.size(); ia++) {
            new_positions[ia].x = positions_[ia].z;
            new_positions[ia].y = positions_[ia].x;
            new_positions[ia].z = positions_[ia].y;
        }
    } else if (new_coord_sequence == "YZX") {
        for (unsigned int ia = 0; ia < positions_.size(); ia++) {
            new_positions[ia].x = positions_[ia].y;
            new_positions[ia].y = positions_[ia].z;
            new_positions[ia].z = positions_[ia].x;
        }
    } else {
        error("Cannot rotate coordinate axes to %s! Possible choices are ZXY and YZX.",
              new_coord_sequence.c_str());
    }
    
    positions_.swap(new_positions);
}


void System::setUnitCell(const vector<Vec3d>& cell_vectors) {
    for (int i = 0; i < 3; i++) {
        cell_matrix_.at(0, i) = cell_vectors[i].x;
        cell_matrix_.at(1, i) = cell_vectors[i].y;
        cell_matrix_.at(2, i) = cell_vectors[i].z;
    }
    tip_pbc_ = true;
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


void System::setDummyXY(double x, double y) {
    real_tip_xy_.x = x;
    real_tip_xy_.y = y;
    
    // If periodic boundary conditions are used, make sure tip position is inside the unit cell
    if (tip_pbc_) {
        Vec3d temp_position = Vec3d(x, y, 0);
        Vec3d cell_norm_position = cell_matrix_.inverse().multiply(temp_position - offset_); // position in basis of unit cell vectors
        cell_norm_position.x -= floor(cell_norm_position.x);
        cell_norm_position.y -= floor(cell_norm_position.y);
        Vec3d pbc_position = cell_matrix_.multiply(cell_norm_position) + offset_;
        
        positions_[0].x = pbc_position.x;
        positions_[0].y = pbc_position.y;
        positions_[1].x = pbc_position.x;
        positions_[1].y = pbc_position.y;
    }
    else {
        positions_[0].x = x;
        positions_[0].y = y;
        positions_[1].x = x;
        positions_[1].y = y;
    }
}
