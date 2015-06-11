#pragma once

#include "simulation.hpp"

int checkForComments(char *line);
int type2num(char *atom);
void parseCommandLine(int argc, char *argv[], Simulation& simulation);
void readInputFile(Simulation& simulation);
void readXYZFile(Simulation& simulation);
void setSystemZ(Simulation& simulation);
void centerSystem(Simulation& simulation);
void readParameterFile(Simulation& simulation);
void parseInteractions(Simulation& simulation, FILE* fp);
