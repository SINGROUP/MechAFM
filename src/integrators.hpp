#pragma once

// Defines the types of integrators available
enum IntegratorType {
    EULER,
    MIDPOINT,
    RK4
};

class System;

// Evolve the system by dt based on Euler integration
void eulerStep(System& system, const double dt);
// Evolve the system by dt based on midpoint integration
void midpointStep(System& system, const double dt);
// Evolve the system by dt based on Runge-Kutta 4 integration
void rk4Step(System& system, const double dt);
