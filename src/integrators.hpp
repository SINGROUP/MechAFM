#pragma once

enum IntegratorType {
    EULER,
    RK4
};

class System;

void eulerStep(System& system);
void rk4Step(System& system);
