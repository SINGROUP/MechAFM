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
    void interpolate(const Vec3d& position, Vec3d& force, double& energy) const;

    Vec3i grid_points_;
    Vec3d spacing_;
    Vec3d offset_;
    const int border_ = 3;  // How many extra samples will we take outside edges of the simulation area
    vector<Vec3d> forces_;
    vector<double> energies_;

 private:
    Vec3i getGridPoint(const Vec3d& pos) const;
    int getGridPointIndex(const Vec3i& grid_point) const;
};
