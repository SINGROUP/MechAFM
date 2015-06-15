#pragma once

#include "globals.hpp"
#include "simulation.hpp"

void error(Simulation& simulation, char *message, ...);
void warning(Simulation& simulation, char *message, ...);
void pretty_print(Simulation& simulation, char *message, ...);
void debugline(Simulation& simulation, int proc, char *message, ...);
void dumpToFiles(Simulation& simulation, BUFFER *sendbuf, BUFFER *recvbuf, int bufsize);
