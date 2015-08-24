#include "minimiser.hpp"

#include <algorithm>
#include <cmath>

#include "integrators.hpp"
#include "interactions.hpp"
#include "messages.hpp"
#include "simulation.hpp"
#include "system.hpp"

using namespace std;

// Checks if the tip force and/or energy have converged within acceptable limits
bool checkConvergence(const Vec3d& tip_force, const double& tip_e_diff,
                                              const InputOptions& options) {
    switch(options.minterm) {
        case MIN_E:
            if (abs(tip_e_diff) < options.etol) {
                return true;
            }
            break;
        case MIN_F:
            if (tip_force.len() < options.ftol) {
                return true;
            }
            break;
        case MIN_EF:
            if (abs(tip_e_diff) < options.etol && tip_force.len() < options.ftol) {
                return true;
            }
            break;
        default:
            error("Invalid minimisation term!");
    }
    return false;
}

int SDMinimisation(System& system, const InputOptions& options) {
    double prev_tip_e = -10e6;
    int n = 1;
    for (; n < options.maxsteps; ++n) {
        // Zero the forces and energies for each step
        if (n != 0) {
            prev_tip_e = system.energies_[1];
        }
        fill(system.forces_.begin(), system.forces_.end(), Vec3d(0));
        fill(system.energies_.begin(), system.energies_.end(), 0);

        // Evaluate all the interactions
        for (const auto& interaction : *system.interactions_) {
            interaction->eval(system.positions_, system.forces_, system.energies_);
        }
        if (checkConvergence(system.forces_[1], system.energies_[1] - prev_tip_e, options)) {
            break;
        }
        for (int i = 0; i < system.n_atoms_; ++i) {
            // Update the atom only if it's not fixed
            if (system.fixed_[i] != 1) {
                system.positions_[i] += options.dt * system.forces_[i];
            }
        }
    }
    return n;
}

int FIREMinimisation(System& system, const InputOptions& options) {
    // Initialize the minimisation variables
    const int n_min = 5;
    const double f_inc = 1.1;
    const double f_dec = 0.5;
    const double f_a = 0.99;
    const double a_start = 0.1;
    const double dt_max = 10 * options.dt;
    double alpha = a_start;
    double dt = options.dt;
    int n_non_neg = 0;
    double prev_tip_e = -10e6;
    int n_tot = 1;
    for (; n_tot < options.maxsteps; ++n_tot) {
        // Zero the forces and energies for each step
        if (n_tot != 0) {
            prev_tip_e = system.energies_[1];
        }
        fill(system.forces_.begin(), system.forces_.end(), Vec3d(0));
        fill(system.energies_.begin(), system.energies_.end(), 0);

        switch (options.integrator_type) {
            case EULER:
                eulerStep(system, dt);
                break;
            case MIDPOINT:
                midpointStep(system, dt);
                break;
            case RK4:
                rk4Step(system, dt);
                break;
            default:
                error("Unimplemented integrator type!");
        }

        if (checkConvergence(system.forces_[1], system.energies_[1] - prev_tip_e, options)) {
            break;
        }

        // Evaluate how we want to change the time step
        double p = 0;
        for (int i = 0; i < system.n_atoms_; ++i) {
            p += system.velocities_[i].dot(system.forces_[i]);
            if (system.fixed_[i] != 1) {
                system.velocities_[i] = (1 - alpha) * system.velocities_[i]
                    + alpha * system.forces_[i].normalized() * system.velocities_[i].len();
            }
        }
        if (p < 0) {
            fill(system.velocities_.begin(), system.velocities_.end(), Vec3d(0));
            n_non_neg = 0;
            dt *= f_dec;
            alpha = a_start;
        } else {
            if (n_non_neg > n_min) {
                dt = min(dt * f_inc, dt_max);
                alpha *= f_a;
            }
            n_non_neg++;
        }
    }
    return n_tot;
}
