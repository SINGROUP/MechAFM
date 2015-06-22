#include "system.hpp"


// Initialize the state vectors and make sure dummy and tip are included
void System::initialize(int n_atoms) {

    n_atoms_ = n_atoms + 2;  // Add tip and dummy to the atom count
    positions_.reserve(n_atoms_);
    velocities_.reserve(n_atoms_);
    forces_.reserve(n_atoms_);
    charges_.reserve(n_atoms_);
    masses_.reserve(n_atoms_);
    fixed_.reserve(n_atoms_);
    types_.reserve(n_atoms_);

    for (int i = 0; i < n_atoms_; ++i) {
        positions_.push_back(Vec3d());
        velocities_.push_back(Vec3d());
        forces_.push_back(Vec3d());
        charges_.push_back(0);
        masses_.push_back(0);
        fixed_.push_back(0);
        types_.push_back("");
    }
}

void System::setDummyXY(double x, double y) {
}

void System::setDummyZ(double z) {
}
