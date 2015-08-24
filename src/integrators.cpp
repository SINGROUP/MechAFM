#include "integrators.hpp"

#include <vector>

#include "system.hpp"
#include "vectors.hpp"

void eulerStep(System& system, const double dt) {
    // Evaluate all the interactions
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(system.positions_, system.forces_, system.energies_);
    }

    for (int i = 0; i < system.n_atoms_; ++i) {
        // Update the atom only if it's not fixed
        if (system.fixed_[i] != 1) {
            system.positions_[i] += dt * system.velocities_[i];
            system.velocities_[i] += dt * system.forces_[i] / system.masses_[i];
        }
    }
}

void midpointStep(System& system, const double dt) {
    // Initialize all the intermediate state vectors
    vector<Vec3d> p2(system.positions_);
    vector<Vec3d> v2(system.velocities_);
    vector<Vec3d> f1(system.n_atoms_, Vec3d(0)), f2(system.n_atoms_, Vec3d(0));
    vector<double> e1(system.n_atoms_, 0), e2(system.n_atoms_, 0);

    // Step 1
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(system.positions_, f1, e1);
    }
    for (int i = 0; i < system.n_atoms_; ++i) {
        if (system.fixed_[i] != 1) {
            p2[i] += dt/2 * system.velocities_[i];
            v2[i] += dt/2 * f1[i] / system.masses_[i];
        }
    }

    // Step 2
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(p2, f2, e2);
    }
    // Update the system
    for (int i = 0; i < system.n_atoms_; ++i) {
        if (system.fixed_[i] != 1) {
            system.positions_[i] += dt * v2[i];
            system.velocities_[i] += dt * f2[i] / system.masses_[i];
        }
        system.forces_[i] = f2[i];
        system.energies_[i] = e2[i];
    }
}

void rk4Step(System& system, const double dt) {
    // Initialize all the intermediate state vectors
    vector<Vec3d> p2(system.positions_), p3(system.positions_), p4(system.positions_);
    vector<Vec3d> v2(system.velocities_), v3(system.velocities_), v4(system.velocities_);
    vector<Vec3d> f1(system.n_atoms_, Vec3d(0)), f2(system.n_atoms_, Vec3d(0));
    vector<Vec3d> f3(system.n_atoms_, Vec3d(0)), f4(system.n_atoms_, Vec3d(0));
    vector<double> e1(system.n_atoms_, 0), e2(system.n_atoms_, 0);
    vector<double> e3(system.n_atoms_, 0), e4(system.n_atoms_, 0);

    // Step 1
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(system.positions_, f1, e1);
    }
    for (int i = 0; i < system.n_atoms_; ++i) {
        if (system.fixed_[i] != 1) {
            p2[i] += dt/2 * system.velocities_[i];
            v2[i] += dt/2 * f1[i] / system.masses_[i];
        }
    }

    // Step 2
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(p2, f2, e2);
    }
    for (int i = 0; i < system.n_atoms_; ++i) {
        if (system.fixed_[i] != 1) {
            p3[i] += dt/2 * v2[i];
            v3[i] += dt/2 * f2[i] / system.masses_[i];
        }
    }

    // Step 3
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(p3, f3, e3);
    }
    for (int i = 0; i < system.n_atoms_; ++i) {
        if (system.fixed_[i] != 1) {
            p4[i] += dt * v3[i];
            v4[i] += dt * f3[i] / system.masses_[i];
        }
    }

    // Step 4
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(p4, f4, e4);
    }
    // Update the system
    for (int i = 0; i < system.n_atoms_; ++i) {
        if (system.fixed_[i] != 1) {
            system.positions_[i] += dt/6 * (system.velocities_[i] + 2*v2[i] + 2*v3[i] + v4[i]);
            system.velocities_[i] += dt/6 * (f1[i] + 2*f2[i] + 2*f3[i] + f4[i]) / system.masses_[i];
        }
        system.forces_[i] = (f1[i] + 2*f2[i] + 2*f3[i] + f4[i]) / 6;
        system.energies_[i] = (e1[i] + 2*e2[i] + 2*e3[i] + e4[i]) / 6;
    }
}
