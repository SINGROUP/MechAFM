#pragma once

enum IntegratorType {
    EULER,
    RK4
};

class System;

void eulerStep(System& system, const double dt);
void rk4Step(System& system, const double dt);
