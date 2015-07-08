#pragma once

#include <memory>
#include <string>
#include <vector>

#include "globals.hpp"
#include "interactions.hpp"
#include "vectors.hpp"

using namespace std;

class System {
 public:
    System(): n_atoms_(0), n_types_(0), interactions_(nullptr) {};
    ~System() {};
    void initialize(int n_atoms);
    OutputData getOutput() const;
    void evalTipSurfaceForces(Vec3d& force, double& energy) const;
    void makeXYZFile() const;
    void setTipDummyDistance(double d) {tip_dummy_d_ = d;}
    void setDummyXY(double x, double y) {
        positions_[0].x = x;
        positions_[0].y = y;
        positions_[1].x = x;
        positions_[1].y = y;
    }
    void setDummyZ(double z) {
        positions_[0].z = z;
        positions_[1].z = z - tip_dummy_d_;
    }

    int n_atoms_;
    int n_types_;
    vector<unique_ptr<Interaction>>* interactions_;

    // Vectors holding the system state
    // index 0 = dummy and index 1 = tip
    vector<Vec3d> positions_;
    vector<Vec3d> velocities_;
    vector<double> charges_;
    vector<double> masses_;
    vector<bool> fixed_;
    vector<string> types_;

 private:
    double tip_dummy_d_;
};
