#pragma once

#include "simulation.hpp"

bool checkForComments(char *line);
void parseCommandLine(int argc, char *argv[], Simulation& simulation);
void readInputFile(Simulation& simulation);
void readXYZFile(Simulation& simulation);
void readParameterFile(Simulation& simulation);
void readFlexibleParameters(Simulation& simulation, FILE* fp);
