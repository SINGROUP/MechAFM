#include "integrators.hpp"

#include "system.hpp"

void eulerStep(System& system, const double dt) {
    // Evaluate all the interactions
    for (const auto& interaction : *system.interactions_) {
        interaction->eval(system.positions_, system.forces_, system.energies_);
    }

    for (int i = 0; i < system.n_atoms_; ++i) {
        // Update the atom only if it's not fixed
        if (!system.fixed_[i]) {
            // printf("%d moved\n", i);
            system.positions_[i] += dt * system.velocities_[i];
            system.velocities_[i] += dt * system.forces_[i]; // / system.masses_[i];
        }
    }
}

void rk4Step(System& system, const double dt) {
    (void)system;
}
