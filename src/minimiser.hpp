#pragma once

#include <memory>
#include <vector>

#include "integrators.hpp"
#include "interactions.hpp"

using namespace std;

enum MinimiserType {
    STEEPEST_DESCENT,
    FIRE
};

struct InputOptions;
class Simulation;
class System;

int SDMinimisation(System& system, const InputOptions& options);
int FIREMinimisation(System& system, const InputOptions& options);
