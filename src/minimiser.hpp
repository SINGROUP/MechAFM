#pragma once

#include <memory>
#include <vector>

#include "system.hpp"

using namespace std;

// Defines the types of minimisers available
enum MinimiserType {
    STEEPEST_DESCENT,
    FIRE
};

struct InputOptions;

// Minimise the system based on criteria given by options with
// Steepest Descent minimisation
int SDMinimisation(System& system, const InputOptions& options);
// Minimise the system based on criteria given by options with
// FIRE minimisation
int FIREMinimisation(System& system, const InputOptions& options);
