#include "force_grid.hpp"

#include <chrono>
#include <cmath>
#include <iostream>

#include "interactions.hpp"
#include "messages.hpp"
#include "simulation.hpp"

using namespace std;


ForceGrid::ForceGrid() {
    is_periodic_ = false;
    is_orthogonal_coord_ = false;
    n_grid_ = Vec3i(0);
    basis_vectors_.assign(3, Vec3d(0));
    offset_ = Vec3d(0);
}


void ForceGrid::setNGrid(const Vec3i& n_grid) {
    n_grid_.x = n_grid.x;
    n_grid_.y = n_grid.y;
    n_grid_.z = n_grid.z;
}


void ForceGrid::setSpacing(const Vec3d& spacing) {
    is_orthogonal_coord_ = true;
    basis_vectors_.assign(3, Vec3d(0));
    basis_vectors_[0].x = spacing.x;
    basis_vectors_[1].y = spacing.y;
    basis_vectors_[2].z = spacing.z;
}


void ForceGrid::setBasisVectors(const vector<Vec3d>& basis_vectors) {
    is_orthogonal_coord_ = false;
    for (int i = 0; i < 3; i++) {
        basis_vectors_[i].x = basis_vectors[i].x;
        basis_vectors_[i].y = basis_vectors[i].y;
        basis_vectors_[i].z = basis_vectors[i].z;
    }
}


void ForceGrid::setOffset(const Vec3d& offset) {
    offset_.x = offset.x;
    offset_.y = offset.y;
    offset_.z = offset.z;
}


void ForceGrid::swapForceValues(vector<Vec3d>& forces_to_swap) {
    forces_.swap(forces_to_swap);
}


void ForceGrid::swapForceValues(DataGrid<Vec3d>& forces_to_swap) {
    forces_to_swap.swapValues(forces_);
}


void ForceGrid::swapEnergyValues(vector<double>& energies_to_swap) {
    energies_.swap(energies_to_swap);
}


void ForceGrid::swapEnergyValues(DataGrid<double>& energies_to_swap) {
    energies_to_swap.swapValues(energies_);
}


void ForceGrid::interpolate(const Vec3d& position, Vec3d& force, double& energy) const {
    // What is the nearest grid point matching position of the tip?
    // Note: grid_point may have negative values if force grid is periodic.
    //       getGridPointIndex() handles periodicity
    Vec3i grid_point = getGridPoint(position);

    // Find the surrounding grid points as array indices and retrieve the force values
    Vec3i current_grid_point;
    vector<Vec3d> force_samples(8);
    vector<double> energy_samples(8);
    for (int i = 0; i < 2; ++i) {
        current_grid_point.x = grid_point.x + i;
        for (int j = 0; j < 2; ++j) {
            current_grid_point.y = grid_point.y + j;
            for (int k = 0; k < 2; ++k) {
                current_grid_point.z = grid_point.z + k;
                // Get the array index
                int index = getGridPointIndex(current_grid_point);
                // Retrieve and store the force
                int sample = i*4 + j*2 + k;
                force_samples[sample] = forces_[index];
                energy_samples[sample] = energies_[index];
            }
        }
    }

    // This is the order of the coefficients in the array
    // 000 - 001 - 010 - 011 - 100 - 101 - 110 - 111

    // The trilinear interpolation below is based on:
    // http://en.wikipedia.org/wiki/Trilinear_interpolation

    // How far are we inside the cube
    Vec3d d;
    d.x = (position.x - offset_.x - grid_point.x * basis_vectors_[1].x) / basis_vectors_[1].x;
    d.y = (position.y - offset_.y - grid_point.y * basis_vectors_[2].y) / basis_vectors_[2].y;
    d.z = (position.z - offset_.z - grid_point.z * basis_vectors_[3].z) / basis_vectors_[3].z;

    // Construct the force
    Vec3d f00, f01, f10, f11, f0, f1;
    f00 = force_samples[0] * (1 - d.x) + force_samples[4] * d.x;
    f01 = force_samples[1] * (1 - d.x) + force_samples[5] * d.x;
    f10 = force_samples[2] * (1 - d.x) + force_samples[6] * d.x;
    f11 = force_samples[3] * (1 - d.x) + force_samples[7] * d.x;
    f0 = f00 * (1 - d.y) + f10 * d.y;
    f1 = f01 * (1 - d.y) + f11 * d.y;
    force = f0 * (1 - d.z) + f1 * (d.z);

    // Construct the energy
    double e00, e01, e10, e11, e0, e1;
    e00 = energy_samples[0] * (1 - d.x) + energy_samples[4] * d.x;
    e01 = energy_samples[1] * (1 - d.x) + energy_samples[5] * d.x;
    e10 = energy_samples[2] * (1 - d.x) + energy_samples[6] * d.x;
    e11 = energy_samples[3] * (1 - d.x) + energy_samples[7] * d.x;
    e0 = e00 * (1 - d.y) + e10 * d.y;
    e1 = e01 * (1 - d.y) + e11 * d.y;
    energy = e0 * (1 - d.z) + e1 * (d.z);
    
}


Vec3i ForceGrid::getGridPoint(const Vec3d& pos) const {
    Vec3i grid_point;
    
    if (is_orthogonal_coord_) {
        grid_point.x = floor((pos.x - offset_.x) / basis_vectors_[0].x);
        grid_point.y = floor((pos.y - offset_.y) / basis_vectors_[1].y);
        grid_point.z = floor((pos.z - offset_.z) / basis_vectors_[2].z);
    }
    else {
        Vec3d pos_in_grid_basis;
    }
    
    // If grid_point is outside the grid inform the user and take the edge point instead.
    // If force grid is periodic, ignore the point being outside the grid.
    if (not is_periodic_) {
        // Check if grid_point.x is outside the grid
        if (grid_point.x < 0) {
            warning("Position outside of grid borders %f, %f, %f!",
                    pos.x, pos.y, pos.z);
            grid_point.x = 0;
        }
        else if (grid_point.x >= n_grid_.x) {
            warning("Position outside of grid borders %f, %f, %f!",
                    pos.x, pos.y, pos.z);
            grid_point.x = n_grid_.x - 1;
        }
        
        // Check if grid_point.y is outside the grid
        if (grid_point.y < 0) {
            warning("Position outside of grid borders %f, %f, %f!",
                    pos.x, pos.y, pos.z);
            grid_point.y = 0;
        }
        else if (grid_point.y >= n_grid_.y) {
            warning("Position outside of grid borders %f, %f, %f!",
                    pos.x, pos.y, pos.z);
            grid_point.y = n_grid_.y - 1;
        }
        
        // Check if grid_point.z is outside the grid
        if (grid_point.z < 0) {
            warning("Position outside of grid borders %f, %f, %f!",
                    pos.x, pos.y, pos.z);
            grid_point.z = 0;
        }
        else if (grid_point.z >= n_grid_.z) {
            warning("Position outside of grid borders %f, %f, %f!",
                    pos.x, pos.y, pos.z);
            grid_point.z = n_grid_.z - 1;
        }
    }
    
    return grid_point;
}


int ForceGrid::getGridPointIndex(const Vec3i& grid_point) const {
    int index = 0;
    
    // Check if the grid_point is outside the grid and react accordingly
    if (grid_point.x < 0 || grid_point.y < 0 || grid_point.z < 0
       || grid_point.x >= n_grid_.x || grid_point.y >= n_grid_.y
       || grid_point.z >= n_grid_.z) {
        
        if (is_periodic_) {
            Vec3i pbc_grid_point = pbcGridPoint(grid_point);
            index = pbc_grid_point.x * n_grid_.y * n_grid_.z
                    + pbc_grid_point.y * n_grid_.z + pbc_grid_point.z;
        }
        else {
            error("Invalid grid point %d, %d, %d (limit: %d, %d, %d)",
                grid_point.x, grid_point.y, grid_point.z,
                n_grid_.x, n_grid_.y, n_grid_.z);
        }
    }
    
    else {
        index = grid_point.x * n_grid_.y * n_grid_.z
                + grid_point.y * n_grid_.z + grid_point.z;
    }
    
    return index;
}


Vec3i ForceGrid::pbcGridPoint(const Vec3i& grid_point) const {
    Vec3i pbc_grid_point = grid_point;
    int n;
    
    // Check if grid_point.x is outside the grid
    if ((grid_point.x < 0) || (grid_point.x >= n_grid_.x)) {
        n = floor(double(grid_point.x)/n_grid_.x);
        pbc_grid_point.x -= n*n_grid_.x;
    }
    
    // Check if grid_point.y is outside the grid
    if ((grid_point.y < 0) || (grid_point.y >= n_grid_.y)) {
        n = floor(double(grid_point.y)/n_grid_.y);
        pbc_grid_point.y -= n*n_grid_.y;
    }
    
    // Check if grid_point.z is outside the grid
    if ((grid_point.z < 0) || (grid_point.z >= n_grid_.z)) {
        n = floor(double(grid_point.z)/n_grid_.z);
        pbc_grid_point.z -= n*n_grid_.z;
    }
    
    return pbc_grid_point;
}
