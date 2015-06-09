#pragma once

#include "simulation.hpp";

int checkForComments(char *line);
int type2num(char *atom);
void parseCommandLine(int argc, char *argv[], Simulation& simulation);
void readInputFile(Simulation& simulation);
void readXYZFile(Simulation& simulation);
void readParameterFile(Simulation& simulation);
