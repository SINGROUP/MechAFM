#include "force_grid.hpp"

#include <cmath>

#include "interactions.hpp"
#include "messages.hpp"
#include "simulation.hpp"

using namespace std;

Vec3i ForceGrid::getGridPoint(const Vec3d& pos) const {
    Vec3i grid_point;
    grid_point.x = floor((pos.x - offset_.x) / spacing_.x);
    grid_point.y = floor((pos.y - offset_.y) / spacing_.y);
    grid_point.z = floor((pos.z - offset_.z) / spacing_.z);
    if (grid_point.x < 0) {
        warning("Position outside of grid borders %f, %f, %f!",
                pos.x, pos.y, pos.z);
        grid_point.x = 0;
    }
    if (grid_point.y < 0) {
        warning("Position outside of grid borders %f, %f, %f!",
                pos.x, pos.y, pos.z);
        grid_point.y = 0;
    }
    if (grid_point.z < 0) {
        warning("Position outside of grid borders %f, %f, %f!",
                pos.x, pos.y, pos.z);
        grid_point.z = 0;
    }
    if (grid_point.x >= grid_points_.x) {
        warning("Position outside of grid borders %f, %f, %f!",
                pos.x, pos.y, pos.z);
        grid_point.x = grid_points_.x - 1;
    }
    if (grid_point.y >= grid_points_.y) {
        warning("Position outside of grid borders %f, %f, %f!",
                pos.x, pos.y, pos.z);
        grid_point.y = grid_points_.y - 1;
    }
    if (grid_point.z >= grid_points_.z) {
        warning("Position outside of grid borders %f, %f, %f!",
                pos.x, pos.y, pos.z);
        grid_point.z = grid_points_.z - 1;
    }
    return grid_point;
}

int ForceGrid::getGridPointIndex(const Vec3i& grid_point) const {
    if (grid_point.x < 0 || grid_point.y < 0 || grid_point.z < 0
       || grid_point.x >= grid_points_.x || grid_point.y >= grid_points_.y
       || grid_point.z >= grid_points_.z) {
        error("Invalid grid point %d, %d, %d (limit: %d, %d, %d)",
                grid_point.x, grid_point.y, grid_point.z,
                grid_points_.x, grid_points_.y, grid_points_.z);
    }
    int index = grid_point.x * grid_points_.y * grid_points_.z
                + grid_point.y * grid_points_.z + grid_point.z;
    return index;
}

void ForceGrid::interpolate(const Vec3d& position, Vec3d& force, double& energy) const {
    // What is the nearest grid point to the current position of the tip?
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

    // How far are we from the corner
    Vec3d d;
    d.x = (position.x - offset_.x - grid_point.x * spacing_.x) / spacing_.x;
    d.y = (position.y - offset_.y - grid_point.y * spacing_.y) / spacing_.y;
    d.z = (position.z - offset_.z - grid_point.z * spacing_.z) / spacing_.z;

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
