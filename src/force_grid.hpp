#pragma once

#include <memory>
#include <vector>

#include "vectors.hpp"

using namespace std;

struct InputOptions;
class Interaction;

class ForceGrid {
 public:
    ForceGrid() {};
    ~ForceGrid() {};
    // Calculates the interpolated force and energy at the given position
    void interpolate(const Vec3d& position, Vec3d& force, double& energy) const;

    Vec3i grid_points_;  // The number of grid points in each dimension
    Vec3d spacing_;  // Distance between grid points in each dimension
    Vec3d offset_;  // The real position of grid point (0, 0, 0)
    const double edge_ = 1.5;  // How wide of an edge do we have around the simulation area
    vector<Vec3d> forces_;  // List of force samples
    vector<double> energies_;  // List of energy samples

 private:
    // Returns the grid point matching for the given position
    Vec3i getGridPoint(const Vec3d& pos) const;
    // Returns the matching array index for the given grid point
    int getGridPointIndex(const Vec3i& grid_point) const;
};
