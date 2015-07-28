#pragma once

#include <memory>
#include <string>
#include <vector>

#include "interactions.hpp"
#include "vectors.hpp"

using namespace std;

class System {
 public:
    System(): n_atoms_(0), interactions_(nullptr) {};
    ~System() {};
    // Initializes the state vectors for given number of surface atoms.
    // Note: n_atoms doesn't include the tip and the dummy!
    void initialize(int n_atoms);
    // Returns the output data for the current state of the system
    OutputData getOutput() const;
    // Evaluates the current force on the tip from surface atoms
    void evalTipSurfaceForces(Vec3d& force, double& energy) const;
    // Writes the current atom positions to a xyz file
    void makeXYZFile(string folder = "") const;
    // Centers the molecule around pos in x, y
    void centerMolecule(const Vec2d pos);
    // Positions the molecule in z to lie on the substrate
    void setMoleculeZ();
    // Sets the tip and dummy initial distance
    void setTipDummyDistance(double d) {tip_dummy_d_ = d;}
    // Returns the tip and dummy initial distance
    double getTipDummyDistance() {return tip_dummy_d_;}
    // Sets dummy x, y coordinates and moves the tip there aswell
    void setDummyXY(double x, double y) {
        positions_[0].x = x;
        positions_[0].y = y;
        positions_[1].x = x;
        positions_[1].y = y;
    }
    // Sets dummy z and moves the tip aswell
    void setDummyZ(double z) {
        positions_[0].z = z;
        positions_[1].z = z - tip_dummy_d_;
    }

    int n_atoms_;  // Count of atoms in the system including the tip and the dummy
    vector<unique_ptr<Interaction>>* interactions_;  // Pointer to the interaction list

    // Vectors holding the system state
    // index 0 = dummy and index 1 = tip
    vector<Vec3d> positions_;
    vector<Vec3d> velocities_;
    vector<Vec3d> forces_;
    vector<double> energies_;
    vector<double> charges_;
    vector<double> masses_;
    vector<int> fixed_;
    vector<string> types_;

 private:
    double tip_dummy_d_;  // Initial distance of the tip and dummy atoms
};
