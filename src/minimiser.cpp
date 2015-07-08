#include "minimiser.hpp"

#include <algorithm>
#include <cmath>

#include "messages.hpp"

using namespace std;

bool checkConvergence(const Vec3d& tip_force, const double& tip_e_diff,
                                              const InputOptions& options) {
    if (options.minterm == MIN_E) {
        if (abs(tip_e_diff) < options.etol) {
            return true;
        }
    } else if (options.minterm == MIN_F) {
        if (tip_force.len() < options.ftol) {
            return true;
        }
    } else if (options.minterm == MIN_EF) {
        if (abs(tip_e_diff) < options.etol && tip_force.len() < options.ftol) {
            return true;
        }
    } else {
        error("Invalid minimisation term!");
    }
    return false;
}

int SDMinimisation(System& system, const InputOptions& options) {
    vector<Vec3d> forces(system.n_atoms_);
    vector<double> energies(system.n_atoms_);
    double prev_tip_e = -10e6;
    int n = 0;
    for (; n < options.maxsteps; ++n) {
        // Zero the forces and energies for each step
        if (n != 0) {
            prev_tip_e = energies[1];
        }
        fill(forces.begin(), forces.end(), Vec3d(0));
        fill(energies.begin(), energies.end(), 0);

        // Evaluate all the interactions
        for (const auto& interaction : *system.interactions_) {
            interaction->eval(system.positions_, forces, energies);
        }
        if (checkConvergence(forces[1], energies[1] - prev_tip_e, options)) {
            break;
        }
        for (int i = 0; i < system.n_atoms_; ++i) {
            // Update the atom only if it's not fixed
            if (!system.fixed_[i]) {
                // printf("%d moved\n", i);
                system.positions_[i] += options.cfac * forces[i];
            }
        }
    }
    return n;
}

int FIREMinimisation(System& system, const InputOptions& options) {
    (void)system;
    (void)options;
    error("FIRE minimiser not yet implemented!");
    return 0;
}
