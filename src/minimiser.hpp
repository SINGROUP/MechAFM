#pragma once

#include <memory>
#include <vector>

#include "system.hpp"

using namespace std;

enum MinimiserType {
    STEEPEST_DESCENT,
    FIRE
};

struct InputOptions;

int SDMinimisation(System& system, const InputOptions& options);
int FIREMinimisation(System& system, const InputOptions& options);
